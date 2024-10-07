st_run = {
 "storages": [
  {

  },
  {

  },
  {

  }
 ],
 
 "run": {

  "center": {
}
}
}



import subprocess
from io import StringIO
import os
import pandas as pd
import time

def check_nas_mount_folder(input_dict):
    unmounted_nas = []
    for item in input_dict.get("storages"):
        mount_location = item.get("mount_location")
        if not os.path.isdir(mount_location):
            unmounted_nas.append(mount_location)
    for nas in unmounted_nas:
        os.makedirs(nas, exist_ok=True)



def is_mounted(folder_path):
    return os.path.ismount(folder_path)


def mount_server(username, server_address, mount_folder, password, version=None):

    if version:
        mount_command = f"sudo mount -t cifs {server_address} {mount_folder} -o username={username},password={password},vers={version}"
    else:
        mount_command = f"sudo mount -t cifs {server_address} {mount_folder} -o username={username},password={password}"

    os.system(mount_command)


def check_reconnect(input_dict):
    
    nas_dict = {}
    mount_locations = []
    for item in input_dict.get("storages"):
        mount_location = item.get("mount_location")
        username = item.get("username")
        password = item.get("password")
        version = item.get("version")
        server_address = item.get("ip_address")
        
        mount_locations.append(mount_location)

        nas_dict[mount_location] = {"username": username, 
                                     "server_address": server_address,
                                     "password": password,
                                     "version": version}


    for folder in mount_locations:

        if not is_mounted(folder):

            print(f"Folder {folder} is not mounted. Reconnecting...")
            username = nas_dict.get(folder).get("username")
            server_address = nas_dict.get(folder).get("server_address")
            version = nas_dict.get(folder).get("version")
            password = nas_dict.get(folder).get("password")

            mount_server(username, server_address, folder, password, version)
            
            time.sleep(5)

            if is_mounted(folder):
                print(f"Successfully connected {server_address} to {folder}")
            else:
                print(f"Failed to connect to {server_address}")
        else:
            print(f"Folder {folder} is already mounted")


def available_nas(input_dict):

    command = "df -h"
    output = subprocess.check_output(command, shell=True, text=True)
    storage = pd.read_csv(StringIO(output), sep='\s+')

    names = []
    orders = []
    locations= []
    mounted = []
    for item in input_dict.get("storages"):
        name = item.get("name")
        order = item.get("order")
        mount_location = item.get("mount_location")
        is_mounted = item.get("ismounted")

        names.append(name)
        orders.append(order)
        locations.append(mount_location)
        mounted.append(is_mounted)


    nas_df = pd.DataFrame({
        "order": orders,
        "name": names,
        "mount_location":locations,
        "is_mounted":mounted 
    })


    nas_df = nas_df.sort_values(by="order", ascending=True)

    mount_locations=[]
    for index, row in nas_df.iterrows():
        ind = storage["Mounted"] == row["mount_location"]
        available_space = storage.loc[ind, "Avail"].values
        if "T" in available_space[0]:
            mount_locations.append(row["mount_location"])
    
    if len(mount_locations) == 0:
    
        return("There is no available space in any of the storage devices")
    
    else:

        return mount_locations[0]


suitable_nas = available_nas(st_run)

suitable_nas

def create_folders_and_path(nas_to_go, input_dict):
    ilyome_path = os.path.join(nas_to_go,"00.Ilyome")
    os.makedirs(ilyome_path, exist_ok=True)       


    run_dict = input_dict.get("run")
    run_id = run_dict.get("id")
    run_name = run_dict.get("run_name")
    center_id = run_dict.get("center").get("id")
    center_name = run_dict.get("center").get("name")

    center_path = os.path.join(ilyome_path, f"{center_id}-{center_name}")

    os.makedirs(center_path, exist_ok=True)

    result_path = os.path.join(center_path,"results")

    os.makedirs(result_path, exist_ok=True)

    runfolder = os.path.join(result_path, f"{run_id}-{run_name}")       
    
    return runfolder


create_folders_and_path(suitable_nas, st_run)



