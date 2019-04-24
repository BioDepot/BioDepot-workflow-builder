# coding: utf-8

from __future__ import absolute_import

from flask import json
from six import BytesIO

from swagger_server.test import BaseTestCase


class TestAgentController(BaseTestCase):
    """AgentController integration test stubs"""

    def test_register_agent(self):
        """Test case for register_agent

        Register Agent
        """
        query_string = [
            ("host_name", "host_name_example"),
            ("cpu_count", 2),
            ("memory", 1),
        ]
        response = self.client.open(
            "/NLPCORE/test2/1.0.0/register-host",
            method="POST",
            query_string=query_string,
        )
        self.assert200(response, "Response body is : " + response.data.decode("utf-8"))


if __name__ == "__main__":
    import unittest

    unittest.main()
