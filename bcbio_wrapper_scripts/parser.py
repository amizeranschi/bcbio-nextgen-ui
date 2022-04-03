# pparse.py

import yaml
import sys
            
def parse_yaml(yml, name=''):
    if isinstance(yml, list):
        for data in yml:
            parse_yaml(data, name)
    elif isinstance(yml, dict):
        if (len(yml) == 1) and not isinstance(yml[list(yml.keys())[0]], list):
            print(str(name+'_'+list(yml.keys())[0]+'='+str(yml[list(yml.keys())[0]]))[1:])
        else:
            for key in yml:
                parse_yaml(yml[key], name+'_'+key)

            
if __name__=="__main__":
    yml = yaml.safe_load(open(sys.argv[1]))
    parse_yaml(yml)