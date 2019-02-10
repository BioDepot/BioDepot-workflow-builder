import json
import sys,os
from threading import Thread

import redis

from .host_registry import HostRegistry


# This agent will find a list of appropriate hosts
# If there are no hosts available then it will return an exception (for now just run locally)
# It will then enqueue the commands onto broker(redis here and the user will provide location and port -
# default will be localhost default Redis port)
# Each host will be sent a run command and the scheduler will return whether the job was successfully started
# redis will handle stdout of each job

def is_redis_available(redis_server, redis_port):
    r = redis.Redis(host=redis_server, port=redis_port)
    try:
        r.get("")  # getting None returns None or throws an exception
    except (redis.exceptions.ConnectionError, redis.exceptions.BusyLoadingError):
        return False
    return True


def enqueue_commands(commands, redis_host, redis_port):
    if not is_redis_available(redis_host, redis_port):
        raise ('server at {}:{} does not exist'.format(redis_host, redis_port))
    # generate queue_id
    queue_id = 'job.{}'.format(os.getpid())
    r = redis.Redis(host=redis_host, port=redis_port)
    for command in commands:
        r.lpush(queue_id, json.dumps(command))
    return queue_id


def send_command(**kwargs):
    HostRegistry.run_command(**kwargs)


def schedule_job(job, cpu_count, memory):  # noqa: E501
    try:
        commands = job['tasks']
    except Exception as e:
        return "Invalid input task, %s" % str(job), 400
    hosts = HostRegistry.get_available_host(core_count=cpu_count, memory=memory)
    if not hosts:
        return "No available agents to process the request.", 400
    else:
        threads = []
        for host in hosts:
            queue_id = enqueue_commands(commands, host.redis_host, host.redis_port)
            t = Thread(target=send_command, kwargs={'host': host, 'queue_id': queue_id, 'redis_host': host.redis_host,
                                                    'redis_port': host.redis_port})
            t.start()
            threads.append(t)
        [t.join() for t in threads]
        [host.free_resources(required_memory=memory, required_cores=cpu_count) for host in hosts]
    # for now return hosts - fr
    return [host.to_object() for host in hosts]
