#!/usr/bin/bash
rsync -av /home/huabei/project/github/RoseTTa-miRP wanghuabei@192.168.0.172:/home/wanghuabei/project --exclude=logs  --exclude=log --exclude=.env --exclude=dataset --exclude=local --exclude=*.tgz --exclude=UniRef30_2020_06 --exclude=RNA --exclude=pdb100_2021Mar03 --exclude=run_RF2NA.sh

