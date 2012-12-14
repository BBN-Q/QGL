class Tek5000():
    def channels(self):
        '''
        The set of empty channels for a Tektronix AWG 5000
        '''
        return {'ch12':{}, 'ch34':{}, 'ch1m1':{}, 'ch1m2':{}, 'ch2m1':{}, 'ch2m2':{}, 'ch3m1':{}, 'ch3m2':{} , 'ch4m1':{}, 'ch4m2':{}}

class Tek7000():
    def channels(self):
        '''
        The set of empty channels for a Tektronix AWG 7000
        '''
        return {'ch12':{}, 'ch34':{}, 'ch1m1':{}, 'ch1m2':{}, 'ch2m1':{}, 'ch2m2':{}, 'ch3m1':{}, 'ch3m2':{} , 'ch4m1':{}, 'ch4m2':{}}


class BBNAPS():
    def channels(self):
        '''
        The set of empty channels for a BBN Arbitrary Pulse Sequencer
        '''
        return {'ch12':{}, 'ch34':{}, 'ch1m1':{}, 'ch2m1':{}, 'ch3m1':{}, 'ch4m1':{}}