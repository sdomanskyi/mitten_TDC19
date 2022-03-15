#!/bin/bash

project_number=
while [ -z $project_number ] #-z means null, so do while project_number is null
do
    echo -n 'What is your Synapse  project number?  '
    read project_number
done

submission_name=
while [ -z $submission_name ]
do
    echo -n 'What is the name of this submission (model name)?  '
    read submission_name
done

tag=
while [ -z $tag ]
do
    echo -n 'What version are you tagging this as?  '
    read tag
done

username=
while [ -z $username ]
do
    echo -n 'What is your Synapse username?  '
    read username
done

password=
while [ -z $password ]
do
    echo -n 'What is your Synapse password?  '
    read password
done



synapse login -u $username  -p $password  --rememberMe
docker build -t docker.synapse.org/syn$project_number/$submission_name:$tag .
docker login docker.synapse.org
docker images
docker push docker.synapse.org/syn$project_number/$submission_name:$tag
