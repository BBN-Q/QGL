import Channels
from Channels import Qubit
from PulsePrimitives import *
from PulseSequencer import Pulse, show, align
import Compiler
from Compiler import compile_to_hardware
from AWG import *

# load the channel params dictionary
Channels.update_channel_info()