import numpy as np
import h5py
import sys
import os.path

def write_to_aps1(fileName, data):
    channelDataFor = np.array([i in data['channelDataFor'] for i in range(1,5)], dtype=np.bool)
    with open(fileName, 'wb') as FID:
        FID.write(b'APS1')                     # target hardware
        FID.write(np.float32(2.2).tobytes())   # Version
        FID.write(channelDataFor.tobytes())    # channelDataFor
        FID.write(np.array(data['miniLLRepeat'], dtype=np.bool).tobytes()) # MiniLLRepeat
        for name in data['channels'].keys():
            FID.write(np.uint8(data['channels'][name]['isIQMode']).tobytes()) # isIQMode
            FID.write(np.uint64(data['channels'][name]['waveformLib'].size).tobytes()) # Length of waveforms
            FID.write(data['channels'][name]['waveformLib'].tobytes()) # Waveforms np.int16
        for name in ['chan_1', 'chan_3']:
            if 'linkListData' in data['channels'][name].keys():
                FID.write(np.uint64(data['channels'][name]['linkListNumKeys']).tobytes()) # numKeys
                FID.write(np.uint64(data['channels'][name]['linkListDataLength']).tobytes()) # numEntries
                for key, dataVec in data['channels'][name]['linkListData'].items():
                    FID.write(key.ljust(32,"#").encode("utf-8")) # Key 32 byte utf-8
                    FID.write(dataVec.tobytes())           

def write_to_aps2(fileName, data):
    instructions = data['instructions']
    wfInfo = {}
    wfInfo[0] = data['chan1']
    wfInfo[1] = data['chan2']

    with open(fileName, 'wb') as FID:
        FID.write(b'APS2')                     # target 
        FID.write(np.float32(data['file_version']).tobytes())   # Version
        FID.write(np.float32(data['fw_version']).tobytes())   # minimum firmware version
        FID.write(np.uint16(2).tobytes())      # number of channels
        # FID.write(np.uint16([1, 2]).tobytes()) # channelDataFor
        FID.write(np.uint64(instructions.size).tobytes()) # instructions length
        FID.write(instructions.tobytes()) # instructions in uint64 form

        #Create the groups and datasets
        for chanct in range(2):
            #Write the waveformLib to file
            if wfInfo[chanct].size == 0:
                #If there are no waveforms, ensure that there is some element
                #so that the waveform group gets written to file.
                #TODO: Fix this in libaps2
                data = np.array([0], dtype=np.int16)
            else:
                data = wfInfo[chanct]
            FID.write(np.uint64(data.size).tobytes()) # waveform data length for channel
            FID.write(data.tobytes())

def get_type(fileName):
    with h5py.File(fileName, 'r') as FID:
        t = FID['/'].attrs['target hardware']
    return t

def read_aps2_from_h5(fileName):
    data = {}
    with h5py.File(fileName, 'r') as FID:
        data["instrument"]   = FID['/'].attrs['target hardware']
        data["file_version"] = FID["/"].attrs["Version"]
        data["fw_version"]   = FID['/'].attrs['minimum firmware version']
        data["chan1"]        = FID['/chan_1/waveforms'].value.flatten()
        data["chan2"]        = FID['/chan_2/waveforms'].value.flatten()
        data["instructions"] = FID['/chan_1/instructions'].value.flatten()
    return data

def read_aps1_from_h5(fileName):
    data = {}
    with h5py.File(fileName, 'r') as FID:
        data["instrument"]     = FID['/'].attrs['target hardware']
        data["file_version"]   = FID["/"].attrs["Version"]
        data["channelDataFor"] = FID["/"].attrs["channelDataFor"]
        data["miniLLRepeat"]   = FID['/'].attrs['miniLLRepeat']
        channels = list(FID['/'].keys())
        data['channels'] = {}
        for channel in channels:
            data['channels'][channel] = {'waveformLib': FID[f'/{channel}/waveformLib'].value.flatten()}
            data['channels'][channel]['isIQMode'] = FID[f'/{channel}'].attrs['isIQMode']
            # print(list(FID[f'/{channel}'].keys()))
            if 'linkListData' in list(FID[f'/{channel}'].keys()):
                data['channels'][channel]['linkListData'] = {}
                for key in FID[f'/{channel}/linkListData'].keys():
                    data['channels'][channel]['linkListNumKeys'] = len(FID[f'/{channel}/linkListData'].keys())
                    data['channels'][channel]['linkListDataLength'] = FID[f'/{channel}/linkListData'].attrs['length']
                    data['channels'][channel]['linkListData'][key] = FID[f'/{channel}/linkListData/{key}'].value.flatten()
    return data

if __name__ == '__main__':
    basename, ext = os.path.splitext(sys.argv[1])
    print(f"Converting {basename+'.h5'}")
    inst = get_type(basename+'.h5')
    if inst == "APS2":
        data = read_aps2_from_h5(basename+'.h5')
        write_to_aps2(basename+".aps2", data)
    if inst == "APS1":
        data = read_aps1_from_h5(basename+'.h5')
        write_to_aps1(basename+".aps1", data)


