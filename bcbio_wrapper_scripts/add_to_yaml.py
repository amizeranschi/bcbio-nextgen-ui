#!/usr/bin/env python
import yaml

document = {}

with open(r'atac-example.yaml') as file:
    documents = yaml.full_load(file)
    document = documents
    document["upload"] = {'dir': '../../final'}

with open(r'./atac-example.yaml', 'w') as file:
    documents = yaml.dump(document, file)
