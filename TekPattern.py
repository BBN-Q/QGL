"""
TekPattern for creating Tek AWG files. Based off of the Matlab version.

Created on Tue May  8 20:39:26 2012

@author: cryan
"""

import struct
import io
import h5py
import numpy as np

MAX_WAVEFORM_VALUE = 2**13-1 #maximum waveform value i.e. 14bit DAC


def write_field(FID, fieldName, data, dataType):
    typeSizes = {'int16':2, 'int32':4, 'double':8, 'uint128':16}
    formatChars = {'int16':'<h', 'int32':'<i', 'double':'<d'}

    if dataType == 'char':
        dataSize = len(data)+1
        data = data+chr(0)
    else:
        dataSize = typeSizes[dataType]
        
    FID.write(struct.pack('<II', len(fieldName)+1, dataSize))
    FID.write(fieldName+chr(0))
    if dataType == 'char':
        FID.write(data)
    elif dataType == 'uint128':
        #struct doesn't support uint128 so write two 64bits
        #there are smarter ways but we really only need this for the fake timestamp
        FID.write(struct.pack('<QQ', 0, data))
    else:
        FID.write(struct.pack(formatChars[dataType], data))
    
def pack_waveform(analog, marker1, marker2):
    '''
    Helper function to convert a floating point analog channel and two logical marker channel to a sequence of 16bit integers.
    % AWG 5000 series binary data format
    % m2 m1 d14 d13 d12 d11 d10 d9 d8 d7 d6 d5 d4 d3 d2 d1
    % 16-bit format with markers occupying left 2 bits followed by the 14 bit
    % analog channel value
    '''

    #Convert decimal shape on [-1,1] to binary on [0,2^14 (16383)] 
    #AWG actually makes 111,111,111,111,10 the 100% output, and
    # 111,111,111,111,11 is one step larger than 100% output so we
    # ignore the one extra positive number and scale from [0,16382]
    analog[analog>1] = 1.0
    analog[analog<-1] = -1.0
    binData = np.uint16( MAX_WAVEFORM_VALUE*analog + MAX_WAVEFORM_VALUE );
    
    binData += 2**14*marker1 + 2**15*marker2
    
    return binData
    
def write_waveform(FID, WFname, WFnumber, data):
    '''
    Helper function to write a waveform
    '''
    numString = str(WFnumber)
    
    write_field(FID, 'WAVEFORM_NAME_'+numString, WFname, 'char')
    
    #Set integer format
    write_field(FID, 'WAVEFORM_TYPE_'+numString, 1, 'int16')

    write_field(FID, 'WAVEFORM_LENGTH_'+numString, data.size, 'int32')
    
    write_field(FID, 'WAVEFORM_TIMESTAMP_'+numString, 0, 'uint128')
    tmpString = 'WAVEFORM_DATA_'+numString+chr(0)
    dataSize = 2*data.size     
    FID.write(struct.pack('<II', len(tmpString), dataSize))
    FID.write(tmpString)
    FID.write(data.tostring())
        

def write_Tek_file(WFs, fileName, seqName, options=None):
    '''
    Main function for writing a AWG format file.
    '''     

    #Set the default options
    #Marker levels default to 1V. 
    if options is None:
        options['markerLevels'] = {}
    for chanct in range(1,5):
        for markerct in range(1,3):
            tmpStr = 'ch{0}m{1}'.format(chanct,markerct)
            if tmpStr not in options['markerLevels']:
                options['markerLevels'][tmpStr] = {}
                options['markerLevels'][tmpStr]['low'] = 0.0
                options['markerLevels'][tmpStr]['high'] = 1.0

    numSeqs = len(WFs['ch1'])
    
    #Open the file
    FID = io.open(fileName, 'wb')

    #Write the necessary MAGIC and VERSION fields
    write_field(FID, 'MAGIC', 5000, 'int16')
    write_field(FID, 'VERSION', 1, 'int16')

    
    #Default to the fastest sampling rate
    write_field(FID, 'SAMPLING_RATE', 1.2e9, 'double')
    
    #Run mode (1 = continuous, 2 = triggered, 3 = gated, 4 = sequence)
    #If we only have one step then there is no sequence
    runMode = 2 if numSeqs == 1 else 4
    write_field(FID, 'RUN_MODE', runMode, 'int16')

    #Default to off state
    write_field(FID, 'RUN_STATE', 0, 'int16')
    
    #Set the reference source (1: internal; 2: external)
    write_field(FID, 'REFERENCE_SOURCE', 2, 'int16')

    #Trigger threshold
    write_field(FID, 'TRIGGER_INPUT_THRESHOLD', 1.0, 'double')    
    
    #Marker's to high/low (1 = amp/offset, 2 = high/low)
    for chanct in range(1,5):
        chanStr = str(chanct)
        write_field(FID, 'CHANNEL_STATE_'+chanStr, 1, 'int16')
        write_field(FID, 'MARKER1_METHOD_'+chanStr, 2, 'int16')
        write_field(FID, 'MARKER1_LOW_'+chanStr, options['markerLevels']['ch'+chanStr+'m1']['low'], 'double')
        write_field(FID, 'MARKER1_HIGH_'+chanStr, options['markerLevels']['ch'+chanStr+'m1']['high'], 'double')
        write_field(FID, 'MARKER2_METHOD_'+chanStr, 2, 'int16')
        write_field(FID, 'MARKER2_LOW_'+chanStr, options['markerLevels']['ch'+chanStr+'m2']['low'], 'double')
        write_field(FID, 'MARKER2_HIGH_'+chanStr, options['markerLevels']['ch'+chanStr+'m2']['high'], 'double')
        
    #If we have only one step then we specify the waveform names
    if numSeqs == 1:
        for chanct in range(1,5):
            write_field(FID, 'OUTPUT_WAVEFORM_NAME_'+str(chanct), seqName+'Ch'+str(chanct)+'001', 'char')
    
    #Now write the waveforms
    for ct in range(numSeqs):
        for chanct in range(1,5):
            chanStr = str(chanct)
            write_waveform(FID, '{0}Ch{1}{2:03d}'.format(seqName, chanct, ct+1), 20+4*ct+chanct, pack_waveform(WFs['ch'+chanStr][ct], WFs['ch'+chanStr+'m1'][ct], WFs['ch'+chanStr+'m2'][ct]))
    
    #Write the sequence table
    for seqct in range(1,numSeqs+1):
        ctStr = str(seqct)
        #We wait for a trigger at every sequence
        write_field(FID, 'SEQUENCE_WAIT_'+ctStr, 1, 'int16')
        write_field(FID, 'SEQUENCE_JUMP_'+ctStr, 0, 'int16')
        write_field(FID, 'SEQUENCE_LOOP_'+ctStr, 1, 'int32')
        
        #If we are on the final one then set the goto back to the beginning
        goto = 1 if seqct == numSeqs else 0
        write_field(FID, 'SEQUENCE_GOTO_'+ctStr, goto, 'int16')
        
        for chanct in range(1,5):
            WFname = '{0}Ch{1}{2:03d}'.format(seqName, chanct, seqct)
            write_field(FID, 'SEQUENCE_WAVEFORM_NAME_CH_'+str(chanct)+'_'+ctStr, WFname, 'char')
            
    FID.close()
    
def read_Tek_file(fileName):
    '''
    Helper function to read in TekAWG h5 dump for plotting.
    '''
    AWGData = {}
    waveformMask = 2**14-1;
    marker1Mask = 2**14;
    marker2Mask = 2**15;
    
    with h5py.File(fileName, 'r') as FID:
        for chanct in range(1,5):
            chanStr = 'ch{0}'.format(chanct)
            marker1Str = 'ch{0}m1'.format(chanct)
            marker2Str = 'ch{0}m2'.format(chanct)
            AWGData[chanStr] = [(1.0/MAX_WAVEFORM_VALUE)*(np.int16(tmpSeq&waveformMask)-MAX_WAVEFORM_VALUE-1) for tmpSeq in FID[chanStr]]
            AWGData[marker1Str] = [tmpSeq&marker1Mask == marker1Mask for tmpSeq in FID[chanStr]];
            AWGData[marker2Str] = [tmpSeq&marker2Mask == marker2Mask for tmpSeq in FID[chanStr]];
    
    return AWGData        
    
if __name__ == '__main__':

    pass
        
    
