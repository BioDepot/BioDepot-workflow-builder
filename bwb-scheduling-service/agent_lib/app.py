import base64
import json
import traceback

from flask import Flask, Response, request


def encode(input_string):
    if isinstance(input_string, str):
        input_string = input_string.encode()
    return base64.b64encode(input_string).decode()


def decode(input_string):
    if isinstance(input_string, str):
        input_string = input_string.encode()
    return base64.b64decode(input_string)


class EndpointAction(object):
    __handler__ = None

    def __init__(self, handler):
        self.__handler__ = handler

    def __call__(self, *args):
        try:
            obj = {'output': encode(self.__handler__())}
            status_code = 200
        except Exception as e:
            obj = {'output': encode(str(e)), 'error': encode(traceback.format_exc())}
            status_code = 500
        return Response(json.dumps(obj), status=status_code)


class FlaskAppWrapper(object):
    app = None
    __hostname__ = None
    __hostport__ = None
    __run_command_fnc__ = None
    __status_fnc__ = None
    __log_fnc__ = None

    def __init__(self, hostname, hostport, run_command_fnc, status_fnc):
        self.__hostname__ = hostname
        self.__hostport__ = hostport
        self.__run_command_fnc__ = run_command_fnc
        self.__status_fnc__ = status_fnc

        app_name = 'agent'
        self.app = Flask(app_name)
        self.add_endpoint(endpoint='/run_command', endpoint_name='run_command', handler=self.run_command,
                          methods=['POST'])
        self.add_endpoint(endpoint='/ping', endpoint_name='ping', handler=ping, methods=['GET'])
        self.add_endpoint(endpoint='/status', endpoint_name='status', handler=self.status, methods=['GET'])


    def run(self):
        self.app.run(host=self.__hostname__, port=self.__hostport__)

    def add_endpoint(self, methods, endpoint=None, endpoint_name=None, handler=None):
        self.app.add_url_rule(endpoint, endpoint_name, EndpointAction(handler), methods=methods)

    def run_command(self):
        data = json.loads(request.data.decode())
        return self.__run_command_fnc__(**data)

    def status(self):
        container_name = request.args.get('container_name')
        if container_name is None:
            raise RuntimeError("get parameter 'container_name' missing.")
        return self.__status_fnc__(container_name=container_name)


def ping():
    return "up"
