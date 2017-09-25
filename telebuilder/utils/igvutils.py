#! /usr/bin/env python

"""
IGV.py

Provides python interface for controlling IGV through a port

"""

import sys
import socket

class IGV:
    def __init__(self, host='0.0.0.0', port=60151, timeout=10, sleep=1, quiet=True):
        self.host = host
        self.port = port
        self.quiet = quiet
        self._connect_socket(timeout, True)
        
        if self.connected:
            self._test()
            if sleep: self.setSleepInterval(int(sleep))
    
    def _connect_socket(self, timeout=10, initial=False):
        if not initial:
            self.socket.close()
        self.socket = socket.socket()
        self.socket.settimeout(timeout)
        try:
            self.socket.connect((self.host, self.port))
            self.connected = True
        except socket.error as e:
            print >>sys.stderr, e
            if e[0] == 61:
                print >>sys.stderr, "Try opening IGV and ensure port is set to %d" % self.port
            self.connected = False
    
    def _test(self):
        if self.send('echo', check_response='echo'):
            if not self.quiet:
                print >>sys.stderr, "Test successful"
        else:
                print >>sys.stderr, "ERROR: Test failed"
    
    def close(self):
        self.socket.close()
    
    def send(self, command, check_response='OK', attempts=5):
        if not self.connected:
            print >>sys.stderr, "ERROR: Not connected to IGV."
            return False
         
        success = False
        attempt = 0
        while success is False:
            attempt += 1
            try:
                self.socket.send('%s\n' % command.strip('\n'))
                response = self.socket.recv(100)
                if response.strip() == check_response:
                    success = True
                else:
                    print >>sys.stderr, 'WARNING: IGV responded "%s"' % response.strip()
            except socket.timeout:
                print >>sys.stderr, "WARNING: socket timed out on attempt %d" % attempt
                self._connect_socket()
            
            if attempt >= attempts:
                break
        
        if not success:
            print >>sys.stderr, "ERROR: Send failed after %d attempts." % attempt
        
        return success
    
    def new(self):
        ''' Create a new session. Unloads all tracks except the default genome annotations
        '''
        self.send("new")
        return self

    def genome(self,genome):
        ''' Selects a genome
        '''
        self.send("genome %s" % genome)
        return self
    
    def load(self, path):
        ''' Loads data or session files. Specify a comma-delimited list of full paths or URLs
        '''
        if not isinstance(path,str):
            self.send('load %s' % ','.join(path))
        else:
            self.send("load %s" % path)
        return self            
    
    def exit(self):
        ''' Exit (close) the IGV application
        '''
        self.send("exit",confirm=False)
    
    def goto(self,locus):
        ''' Scrolls to a single locus or space-delimited list of loci. If a list is provided,
            these loci will be displayed in a split screen view. Use any syntax that is valid
            in the IGV search box.
        '''
        self.send("goto %s" % locus)
    
    def region(self, chr, start, end):
        ''' Not implemented '''
        pass
  
    def maxPanelHeight(self, height):
        ''' Sets the number of vertical pixels (height) of each panel to include in image
        '''
        self.send("maxPanelHeight %d" % height)
          
    def snapshot(self, filename=None):
        ''' Saves a snapshot of the IGV window to an image file. If filename is omitted,
            writes a PNG file with a filename generated based on the locus. If filename is
            specified, the filename extension determines the image file format, which must be
            .png, .jpg, or .svg
        '''
        extensions = set(['png','svg','jpg'])
        if filename is None: self.send('snapshot')
        else:
            assert filename.split('.')[-1] in extensions , 'ERROR: Filename "%s" is invalid' % filename
            self.send('snapshot %s' % filename)
  
    def snapshotDirectory(self,dname='.'):
        ''' Sets the directory in which to write images
        '''
        import os
        assert os.path.isdir(dname)
        self.send('snapshotDirectory %s' % os.path.abspath(dname))
    
    def viewaspairs(self, trackName=None):
        ''' Set the display mode for an alignment track to "View as pairs". trackName is
            optional.
        '''
        if trackName is None: self.send('viewaspairs')
        else: self.send('viewaspairs %s' % trackName)
  
    def squish(self, trackName=None):
        ''' Squish a given trackName. trackName is optional, and if it is not supplied all
            annotation tracks are squished.
        '''
        if trackName is None: self.send('squish')
        else: self.send('squish %s' % trackName)
  
    def collapse(self, trackName=None):
        ''' Collapse a given trackName. trackName is optional, and if it is not supplied all
            annotation tracks are collapsed.
        '''
        if trackName is None: self.send('collapse')
        else: self.send('collapse %s' % trackName)
    
    def expand(self, trackName=None):
        ''' Expand a given trackName. trackName is optional, and if it is not supplied all
            annotation tracks are expanded.
        '''
        if trackName is None: self.send('expand')
        else: self.send('expand %s' % trackName)
        return self        
    
    def sort(self, option, locus):
        ''' Not implemented yet '''
        pass
        
    def preference(self, trackName=None):
        ''' Not implemented '''
        pass
    
    def setSleepInterval(self, ms=1000):
        ''' Sets a delay (sleep) time in milliseconds.  The sleep interval is invoked 
            between successive commands.
        '''
        if self.send('setSleepInterval %d' % int(ms)):
            if not self.quiet:
                print >>sys.stderr, 'Setting sleep interval to %d' % int(ms)


def igv_init(genome, tracks=[]):
    igv = IGV()
    if not igv.connected:
        return None

    igv.new().genome(genome)
    for t in tracks:
        igv.load(t)
    return igv
