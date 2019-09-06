import json
import os,sys

import threading
import redis

from .host_registry import HostRegistry, JobUnit

# This agent will find a list of appropriate hosts
# If there are no hosts available then it will return an exception (for now just run locally)
# It will then enqueue the commands onto broker(redis here and the user will provide location and port -
# default will be localhost default Redis port)
# Each host will be sent a run command and the scheduler will return whether the job was successfully started
# redis will keep track of the list containers for each job and the pids to stop

def is_redis_available(redis_server, redis_port):
    r = redis.Redis(host=redis_server, port=redis_port)
    try:
        r.get("")  # getting None returns None or throws an exception
    except (redis.exceptions.ConnectionError, redis.exceptions.BusyLoadingError):
        return False
    return True


def enqueue_commands(commands, redis_host, redis_port, job_id):
    if not is_redis_available(redis_host, redis_port):
        raise ('server at {}:{} does not exist'.format(redis_host, redis_port))
    try:
        job_queue_id = 'jobs'
        command_queue_id = 'commands.{}'.format(job_id)
        r = redis.Redis(host=redis_host, port=redis_port)
        r.sadd(job_queue_id,str(job_id))
        for command in commands:
            r.lpush(command_queue_id, json.dumps(command))
        return command_queue_id
    except Exception as e:
        r.srem(job_ids,job_id)
        return None
        
def cleanup_redis (job_id,redis_host, redis_port):
    if not is_redis_available(redis_host, redis_port):
        raise ('server at {}:{} does not exist'.format(redis_host, redis_port))
    try:
        job_queue_id='jobs'
        r.srem(job_queue_id,job_id)
        command_queue_id = 'commands.{}'.format(job_id)
        q=r.rpop(command_queue_id) 
        while q :
            q=r.rpop(command_queue_id) 
    except Exception as e:
        r.srem(jobs_ids,job_id)
        return None    

def schedule_job(job_id, job, cpu_count, memory):
    try:
        commands = job['tasks']
        maxWorkers = job ['maxWorkers']
    except Exception as e:
        return "Invalid input task, %s" % str(e), 400
    job_unit = HostRegistry.get_available_host(core_count=cpu_count, memory=memory, maxWorkers=maxWorkers, job_id=job_id)
    if not job_unit.hosts:
        return "No available agents to process the request.", 400
    else:
        redis_host=job_unit.hosts[0].redis_host
        redis_port=job_unit.hosts[0].redis_port
        job_id=job_unit.job_id
        command_queue_id = enqueue_commands(commands,redis_host, redis_port, job_id)
        threads=[]
        if command_queue_id:
            for i in range(len(job_unit.hosts)):
                thread=threading.Thread(target=HostRegistry.run_command,args=(job_unit.hosts[i],command_queue_id,redis_host,redis_port,job_id,str(i+1)))
                thread.daemon = True
                threads.append(thread)
                thread.start()
        for thread in threads:
            thread.join()
        #change this to use redis for bookkeeping so that threads can signal when they are finished
        for host in job_unit.hosts:
            host.free_resources(required_memory=memory, required_cores=cpu_count)
        return 0

def status(job_id):
    try:
        job_id_obj = __ACTIVE_NAMESPACES__[job_id]
    except KeyError:
        return "Namespace does not exist", 404
    status_codes = {}
    for active_hosts in job_id_obj:
        containers = active_hosts['containers']
        host = active_hosts['host']
        for container_name in containers:
            status_code = HostRegistry.status(host, container_name=container_name)
            status_codes[container_name] = status_code
    return status_codes


def log(job_id,worker_id=0):
    return HostRegistry.log(job_id,worker_id)
