import json
from typing import Dict


def read_json(path: str) -> Dict:
    """Read the json file.
    Args
    :path the path to the file
    :return a dictionary of the json file
    """
    if path.endswith(".json"):
        with open(path, "rb") as f:
            content = json.loads(f.read())
            return content
    else:
        print("Not a json file")
        return {}


def save_json(content: Dict, path: str):
    """Save the dictionary into a json file.
    Args
    :param content: the object to save
    :param path: the path where to save. Could have .json or not in the path
    """
    assert type(content) is dict
    if path.endswith(".json"):
        path_to_save = path
    else:
        path_to_save = path + ".json"
    with open(path_to_save, "w") as file:
        json.dump(content, file)
