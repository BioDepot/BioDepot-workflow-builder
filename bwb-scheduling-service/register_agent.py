import argparse
from sys import argv

from agent_lib.agents import DockerLocalAgent

AGENT_TYPES = {
    'local-agent': DockerLocalAgent
}


def main():
    parser = argparse.ArgumentParser(description="Agent utility")
    parser.add_argument("agent_type", type=str, choices=['local-agent'])
    parser.add_argument("--api_url", type=str, required=True)
    parser.add_argument("--hostname", type=str, default='localhost')
    parser.add_argument("--hostport", type=int, default=8090)
    parser.add_argument("--cpu_count", type=int, default=0)
    parser.add_argument("--mem_limit", type=int, default=0)
    parser.add_argument("agent_args", nargs=argparse.REMAINDER)
    parsed_args = parser.parse_args(argv[1:])
    agent_type = AGENT_TYPES[parsed_args.agent_type]

    agent = agent_type(argv=parsed_args.agent_args, api_url=parsed_args.api_url, hostname=parsed_args.hostname,
                       hostport=parsed_args.hostport, cpu_count=parsed_args.cpu_count, mem_limit=parsed_args.mem_limit)
    agent.start_agent()
    agent.register_agent()
    agent.join()


if __name__ == "__main__":
    main()
