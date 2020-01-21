#!/bin/bash

oldlogs=$(find logs -ctime +100 -name "*.log")
tar -uzf logs/oldlogs.tar.gz $oldlogs &&
rm $oldlogs
