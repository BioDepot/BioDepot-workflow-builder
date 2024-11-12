#!/bin/bash
if xmessage -nearmouse -buttons No:1,Yes:0 "Are you sure you want to quit?"; then
    pkill supervisord
else
  xmessage -nearmouse "Quit aborted"
fi