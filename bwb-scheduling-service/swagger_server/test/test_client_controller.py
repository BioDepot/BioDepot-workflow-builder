# coding: utf-8

from __future__ import absolute_import

from flask import json
from six import BytesIO

from swagger_server.models.inventory_item import InventoryItem  # noqa: E501
from swagger_server.test import BaseTestCase


class TestClientController(BaseTestCase):
    """ClientController integration test stubs"""

    def test_schedule_job(self):
        """Test case for schedule_job

        Schedule Job
        """
        query_string = [("command", "command_example"), ("cpu_count", 2), ("memory", 1)]
        response = self.client.open(
            "/NLPCORE/test2/1.0.0/schedule-job", method="GET", query_string=query_string
        )
        self.assert200(response, "Response body is : " + response.data.decode("utf-8"))


if __name__ == "__main__":
    import unittest

    unittest.main()
