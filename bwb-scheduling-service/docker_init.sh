#!/bin/bash

command=$1
args="${@:2}"

if [ "$command" = "service" ]; then
    python3 -m swagger_server $args
elif [ "$command" = "agent" ]; then
    python3 register_agent.py $args
else
    echo "Usage service|agent <args>"
fi

