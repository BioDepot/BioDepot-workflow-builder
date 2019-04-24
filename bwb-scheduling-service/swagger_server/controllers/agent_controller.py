from .host_registry import HostRegistry


def register_agent(
    host_name, host_port, core_count, memory, redis_host, redis_port
):  # noqa: E501
    HostRegistry.register_resource(
        host_name=host_name,
        host_port=host_port,
        core_count=core_count,
        memory=memory,
        redis_host=redis_host,
        redis_port=redis_port,
    )
