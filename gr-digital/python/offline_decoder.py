import math
from gnuradio import gr
import digital_swig
import ofdm_packet_utils
from ofdm_receiver import ofdm_receiver
import gnuradio.gr.gr_threading as _threading
import psk, qam


class ofdm_decode(gr.hier_block2):

    def __init__(self, options, callback=None):

	gr.hier_block2.__init__(self, "ofdm_decode",
				gr.io_signature(1, 1, gr.sizeof_gr_complex), # Input signature
				gr.io_signature(1, 1, gr.sizeof_gr_complex)) # Output signature

        self._rcvd_pktq = gr.msg_queue()          # holds packets from the PHY
        self._modulation = options.modulation
        self._fft_length = options.fft_length
        self._occupied_tones = options.occupied_tones
        self._cp_length = options.cp_length
        self._snr = options.snr
        self._threshold = options.threshold
        self._ra = options.rate_adapt
        print "ra=",self._ra
        print "threshold,",self._threshold

        mods = {"bpsk": 2, "qpsk": 4, "8psk": 8, "qam8": 8, "qam16": 16, "qam64": 64, "qam256": 256}
        arity = mods[self._modulation]

        rot = 1
        if self._modulation == "qpsk":
            rot = (0.707+0.707j)


        # FIXME: pass the constellation objects instead of just the points
        if(self._modulation.find("psk") >= 0):
            constel = psk.psk_constellation(arity)
            rotated_const = map(lambda pt: pt * rot, constel.points())
        elif(self._modulation.find("qam") >= 0):
            constel = qam.qam_constellation(arity)
            rotated_const = map(lambda pt: pt * rot, constel.points())
 
        phgain = 0.25
        frgain = phgain*phgain / 4.0

        self.ofdm_demod = digital_swig.ofdm_partialcomb_receiver(rotated_const, range(arity),
                                                       self._rcvd_pktq,
                                                       self._occupied_tones,
                                                       self._ra,
                                                       phgain, frgain)

        self.connect(gr.file_source(gr.sizeof_gr_complex*self._occupied_tones, "ofdm_frame_sink0_c.dat"), (self.ofdm_demod,0)) 
        self.connect(gr.file_source(gr.sizeof_char, "ofdm_frame_sink1_char.dat"), (self.ofdm_demod,1)) 


    def add_options(normal, expert):
        """
        Adds OFDM-specific options to the Options Parser
        """
        normal.add_option("-m", "--modulation", type="string", default="bpsk",
                          help="set modulation type (bpsk or qpsk) [default=%default]")
        expert.add_option("", "--fft-length", type="intx", default=512,
                          help="set the number of FFT bins [default=%default]")
        expert.add_option("", "--occupied-tones", type="intx", default=200,
                          help="set the number of occupied FFT bins [default=%default]")
        expert.add_option("", "--cp-length", type="intx", default=128,
                          help="set the number of bits in the cyclic prefix [default=%default]")
        expert.add_option("", "--snr", type="float", default=30.0,
                          help="SNR estimate [default=%default]")
        expert.add_option("", "--threshold", type="float", default=1.0,
                          help="cross correlation threshold [default=%default]")
        expert.add_option("", "--rate-adapt", type="intx", default=0,
                          help="rate adaptation scheme [default=%default]")


    # Make a static method to call before instantiation
    add_options = staticmethod(add_options)

    def _print_verbage(self):
        """
        Prints information about the OFDM demodulator
        """
        print "\nOFDM Demodulator:"
        print "Modulation Type: %s"    % (self._modulation)
        print "FFT length:      %3d"   % (self._fft_length)
        print "Occupied Tones:  %3d"   % (self._occupied_tones)
        print "CP length:       %3d"   % (self._cp_length)



class _queue_watcher_thread(_threading.Thread):
    def __init__(self, rcvd_pktq, callback):
        _threading.Thread.__init__(self)
        self.setDaemon(1)
        self.rcvd_pktq = rcvd_pktq
        self.callback = callback
        self.keep_running = True
        self.start()


    def run(self):
        while self.keep_running:
            msg = self.rcvd_pktq.delete_head()
            ok, payload = ofdm_packet_utils.unmake_packet(msg.to_string())
            if self.callback:
                self.callback(ok, payload)


