from gnuradio import gr, blks2
from gnuradio import eng_notation
from gnuradio.eng_option import eng_option
from optparse import OptionParser

from gnuradio import digital

# from current dir
from offline_path import offline_path
from uhd_interface import uhd_receiver
import struct, sys, os
print os.getpid()
raw_input("Press enter to continue")



class my_top_block(gr.top_block):
    def __init__(self, callback, fwd_callback, options):
        gr.top_block.__init__(self)

        if(options.rx_freq is not None):
            self.source = uhd_receiver(options.args,
                                       options.bandwidth,
                                       options.rx_freq, options.rx_gain,
                                       options.spec, options.antenna,
                                       options.verbose)
        elif(options.from_file is not None):
            self.source = gr.file_source(gr.sizeof_gr_complex, options.from_file)
        else:
            self.source = gr.null_source(gr.sizeof_gr_complex)

        # Set up receive path
        # do this after for any adjustments to the options that may
        # occur in the sinks (specifically the UHD sink)
        self.rxpath = offline_path(callback, fwd_callback, options)

        self.connect(self.source, self.rxpath)
	#self.connect(self.rxpath)
 
def main():

    def rx_callback():
            print "rx_callback (wrapper) invoked!"


    def fwd_callback():
            print "fwd_callback (wrapper) invoked!"


    parser = OptionParser(option_class=eng_option, conflict_handler="resolve")
    expert_grp = parser.add_option_group("Expert")
    parser.add_option("","--discontinuous", action="store_true", default=False,
                      help="enable discontinuous")
    parser.add_option("","--from-file", default=None,
                      help="input file of samples to demod")

    parser.add_option("", "--snr", type="eng_float", default=30, help="set the SNR of the channel in dB [default=%default]")

    offline_path.add_options(parser, expert_grp)
    uhd_receiver.add_options(parser)
    digital.ofdm_decode.add_options(parser, expert_grp)

    (options, args) = parser.parse_args ()

    if options.from_file is None:
        if options.rx_freq is None:
            sys.stderr.write("You must specify -f FREQ or --freq FREQ\n")
            parser.print_help(sys.stderr)
            sys.exit(1)


    # build the graph
    tb = my_top_block(rx_callback, fwd_callback, options)

    r = gr.enable_realtime_scheduling()
    if r != gr.RT_OK:
        print "Warning: failed to enable realtime scheduling"

    tb.start()                      # start flow graph
    tb.wait()                       # wait for it to finish


if __name__ == '__main__':
    try:
        main()
    except KeyboardInterrupt:
        pass
