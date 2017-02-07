import struct
import os
import math

PP_start_header='%start_PP_Header%'
PP_start_string='%START_PP_DATA%'
MER_start_string='%START_DATA%'

def read_MER_simulations():
    input = []
    output = []

    #load all the simulations and create the input and output data for cnn
    for ppsim in [x for x in os.listdir('.') if not x.rfind('ppsim')]:
        with open(ppsim,'rb') as fp:
            file_in = fp.read()
        blocks = file_in.split(PP_start_string)

        #Parse the headers for the distribution parameters
        headers = blocks[0].split(PP_start_header)
        for tokens in headers[0].split(';'):
            token = tokens.split(':')

            #there are lots more there but I generated
            # these recordings so I'll just deal with 
            # what I want for now
            if token[0] == 'distribution_parameters':
                parameters = [float(x) for x in token[1][1:-1].split(',')][1:]
            elif token[0] == 'fs':
                sample_rate = int(token[1])
        for tokens in headers[1].split(';'):
            token = tokens.split(':')

            #Same as above
            if token[0] == 'pp_length':
                no_samples = int(token[1])

        #load only the MER simulation I don't need anything else
        data_blocks = blocks[1].split(MER_start_string)
        mers = struct.unpack('d'*no_samples, data_blocks[1])

        #chunk the data into 1s recordings skipping the first 10s to 
        # allow model run in (mainly for stretched exponential). 
        n = 24000
        parameters[0] = math.log(parameters[0])
        parameters[1] = parameters[1] * sample_rate
        for mer in [mers[i:i + n] for i in range(no_samples//2, no_samples - 2*n, n)]:
            input.append(mer)
            output.append(parameters)
    return (input, output)