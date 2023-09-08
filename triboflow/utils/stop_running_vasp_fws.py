import os

from fireworks import LaunchPad

lp_dict = {
    "authsource": "admin",
    "host": "mongodb+srv://surfgen.d3fwalo.mongodb.net",
    "logdir": None,
    "mongoclient_kwargs": {},
    "name": "fireworks_release_test",
    "password": "xEigc3v4GJ52B6XJ",
    "port": 27017,
    "strm_lvl": "INFO",
    "uri_mode": False,
    "user_indices": [],
    "username": "seadmin",
    "wf_user_indices": [],
}

lpad = LaunchPad.from_dict(lp_dict)

wf_id = 1536

wflow = lpad.get_wf_by_fw_id(wf_id)

for fw_id, state in wflow.fw_states.items():
    if state == "RUNNING":
        print(f"Placing STOPCAR {fw_id}...")
        ldir = lpad.get_launchdir(fw_id)
        # check if vasp is running in this directory by checking for an OUTCAR file
        # if not, skip this directory
        if not os.path.exists(ldir + "/OUTCAR"):
            continue
        # if vasp is running, then
        # write a file named STOPCAR in the launch directory. It contains
        # the line "LABORT = .TRUE." and nothing else.
        with open(ldir + "/STOPCAR", "w") as f:
            f.write("LABORT = .TRUE.")
