import re


with open("biodepot/setup.py", "r") as fp:
    matches=re.findall('name="(.*)"',fp.read())

for match in matches:
    print(match)
    
