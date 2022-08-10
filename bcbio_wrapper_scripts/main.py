import eel
import os
import json
import ruamel.yaml
import time
import yaml_to_table

from tkinter import *
from tkinter.filedialog import askopenfilename

eel.init("web")

conf_path = os.getcwd() + '/result.yaml'
report_context_path = os.getcwd() + '/web/data/report_context.json'

# Exposing the python function to javascript

@eel.expose
def save_configuration(config):
    try:
        dict = {}
        for k, v in config.items():
            if v['type'] == 'text' or v['type'] == 'checkbox':
                dict[k] = v['value']
            if v['type'] == 'number':
                dict[k] = int(v['value'])

        with open(conf_path, 'w') as yaml_file:
            yaml = ruamel.yaml.YAML()
            yaml.indent(sequence=4, offset=2)
            yaml.dump(dict, yaml_file)
        return conf_path
    except:
        return ''


@eel.expose
def generate_preview():
    result = yaml_to_table.generatePreviewFromYaml(conf_path)
    return result


@eel.expose
def load_config_file():
    root = Tk()
    root.withdraw()
    root.call('wm', 'attributes', '.', '-topmost', True)
    file_path = askopenfilename()
    root.update()
    if not file_path:
        return 'cancel'
    if file_path.endswith('yml') or file_path.endswith('yaml'):
        with open(file_path, 'r') as stream:
            try:
                yaml = ruamel.yaml.YAML()
                yaml.indent(sequence=4, offset=2)
                data = dict(yaml.load(stream))

                with open(conf_path, 'w') as yaml_file:
                    yaml = ruamel.yaml.YAML()
                    yaml.indent(sequence=4, offset=2)
                    yaml.dump(data, yaml_file)
                return dict(data)
            except:
                return None
    else:
        return None

@eel.expose
def get_report_context_file():
	# Opening JSON file
	f = open(report_context_path)
	
	# returns JSON object as 
	# a dictionary
	data = json.load(f)
	return data

@eel.expose
def run_analysis():
	print("Trigger All scripts")
	try:
		os.system('bash ' + os.getcwd()+ '/deploy.sh ' + conf_path)
		os.system('cp -r ' + os.getcwd()+ '/downstreamAnalysisBulk-RNA-seq ' + os.getcwd() + '/web/images/')
		# time.sleep(10)
		return True
	except:
		return False

	# Return True for success and False for errors

# Start the index.html file
eel.start("index.html")
