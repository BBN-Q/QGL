#!/usr/bin/env python

import sys
import socket
import struct
try:
    from QGL.drivers import APS2Pattern
    haveQGL = True
except:
    haveQGL = False
import numpy as np
from ansicolor import * # from ansicolors

ip = sys.argv[1]
port = 0xbb50

short_time = True

raw_mode = '-r' in sys.argv
tdm_mode = '-tdm' in sys.argv

print("Connecting to APS Debug Port at {0}:{1}".format(ip,port))

PACKET_SIZE=168//8+1

if short_time:
    PACKET_SIZE -= 2

MODE_UNKNOWN = 0
MODE_SEQUENCER = 1
MODE_RAM = 2
MODE_TRIGGER = 3

def bittest(v, b):
    if ((v & (1<<b)) == (1<<b)):
        return red('T')
    else:
        return 'F'

def get_bytes(a, s, e):
    e += s
    return a[s:e], e

def formatHaltBits(haltBits):
    haltSync = bittest(haltBits,5)
    haltCustom = bittest(haltBits,4)
    haltValid = bittest(haltBits,3)
    haltWait = bittest(haltBits,2)
    haltStr = "S:{} V:{} C:{} W:{}".format(haltSync, haltValid, haltCustom, haltWait)
    return haltStr



s = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
s.connect((ip,port))
print("Connected")
print(s)
extraData = b''
while 1:
    data = extraData + s.recv(10*PACKET_SIZE)
    num_packets = len(data) // PACKET_SIZE
    num_excess = len(data) % PACKET_SIZE
    block_size = num_packets*PACKET_SIZE
    extraData = b''
    if num_excess > 0:
        extraData = data[block_size:]
        data = data[:block_size]

    for cnt in range(num_packets):
        start = cnt * PACKET_SIZE
        end = (cnt+1)*PACKET_SIZE
        packet = data[start:end]
        packet_bytes = bytearray(packet)
        if raw_mode:
            print(packet_bytes)
            continue

        if packet[0] ==  0x1:
            mode = MODE_SEQUENCER
        elif packet[0] in [0x2,0x3, 0x4, 0x5]:
            mode = MODE_RAM
        elif packet[0] == 0x6:
            mode = MODE_TRIGGER
        else:
            mode = MODE_UNKNOWN


        if short_time:
            uptime_seconds = packet_bytes[1:3]
            uptime_nanoseconds = packet_bytes[3:7]

            uptime_seconds = struct.unpack(">H", uptime_seconds)[0]
        else:
            uptime_seconds = packet_bytes[1:5]
            uptime_nanoseconds = packet_bytes[5:9]

            uptime_seconds = struct.unpack(">I", uptime_seconds)[0]

        uptime_nanoseconds = struct.unpack(">I", uptime_nanoseconds)[0]

        #print(uptime_seconds, uptime_nanoseconds)

        uptime_nanoseconds = uptime_nanoseconds/1e9
        uptime = uptime_seconds + uptime_nanoseconds

        if short_time:
            start = 8
        else:
            start = 10

        if mode == MODE_SEQUENCER:

            if short_time:
                triggerWord = packet[7]
            else:
                triggerWord = packet[9]

            haltBits = packet[start]

            haltStr = formatHaltBits(haltBits)

            instructionAddr, start = get_bytes(packet, start, 4)
            instructionAddr = bytearray(instructionAddr)
            # clear halt bits, only the first two bits are part of the Addr
            instructionAddr[0] = instructionAddr[0] & 0x3
            instructionAddr = struct.unpack(">I", instructionAddr)[0]

            seq_debug_data, start = get_bytes(packet, start, 8)
            seq_debug_data = np.fromstring(seq_debug_data[::-1], dtype=np.uint64)

            if haveQGL:
                instruction = APS2Pattern.Instruction.unflatten(seq_debug_data, decode_as_tdm = tdm_mode)
            else:
                instruction = ''

            h = "{:016x}".format(seq_debug_data[0]).upper()

            print(red("{:10}".format("Sequencer")), \
                  white(": {:4}".format(instructionAddr)),\
                  blue("0x{}".format(h)), \
                  yellow(" {:6.8f}".format(uptime)), \
                  haltStr, \
                  "T:0x{:x}".format(triggerWord), \
                  " {} ".format(instruction))
            # if instructionAddr > 100:
            #    sys.exit()
        elif mode == MODE_RAM:

            zeros, start = get_bytes(packet, start, 4)
            vram_addr, start = get_bytes(packet, start, 4)
            vram_data, start = get_bytes(packet, start, 4)

            vram_header = struct.unpack(">I", zeros)[0]
            vram_addr = struct.unpack(">I", vram_addr)[0]
            vram_data = struct.unpack(">I", vram_data)[0]

            if packet[0] == 0x2:
                ram_mode = "Valid"

            if packet[0] == 0x3:
                ram_mode = "Write"

            if packet[0] == 0x4:
                ram_mode = "Send"

            if packet[0] == 0x5:
                ram_mode = "Recv"

            print(green("{:10}".format("RAM")), \
                  white(": {:>23}".format(ram_mode)),\
                  yellow(" {:6.8f}".format(uptime)), \
                  "addr = 0x{:08X} data = 0x{:08X}".format(vram_addr, vram_data))

        elif mode == MODE_TRIGGER:

            zeros, start = get_bytes(packet, start, 8)
            syncBits, start = get_bytes(packet,start,1)
            haltBits, start = get_bytes(packet, start, 1)
            triggers, start = get_bytes(packet, start, 1)
            triggerWord, start = get_bytes(packet,start,1)

            for z in zeros:
                if z != 0:
                    print("Error expected 0 not", z)

            haltBits = haltBits[0]
            syncBits = syncBits[0]
            triggers = triggers[0]
            triggerWord = triggerWord[0]

            haltStr = formatHaltBits(haltBits)

            syncWF = (syncBits >> 3) & 0xF
            syncMarker = (syncBits >> 1) & 0x3
            syncMod = syncBits & 0x1


            halt = bittest(syncBits,4)
            pc_jump = bittest(syncBits, 5)
            tready = bittest(syncBits, 6)
            tvalid = bittest(syncBits, 7)

            syncStr = "TV: {} TR: {} PJ: {} H: {} SWF: 0x{} SMa: 0x{} SMo: 0x{}".format(tvalid, tready, pc_jump, halt, syncWF, syncMarker, syncMod)

            trigger = bittest(triggers,1)
            triggerWordValid = bittest(triggers,0)

            triggerStr = "t = {} tWV = {} tW = 0x{:0X}".format(trigger,triggerWordValid, triggerWord)

            print(green("{:10}".format("Trigger")), \
                  white(": {:>23}".format('')),\
                  yellow(" {:6.8f}".format(uptime)), \
                  haltStr, \
                  syncStr, \
                  triggerStr)
        elif mode == MODE_UNKNOWN:
            print(green("{:10}".format("Unknown")),packet_bytes)

s.close()
