#! /usr/bin/env python

"""
IGV.py

Provides python interface for controlling IGV through a port

"""

import sys
import socket

class IGV:
    def __init__(self, host='127.0.0.1', port=60151, timeout=10, sleep=1, quiet=True):
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
        self.socket = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
        # self.socket.settimeout(timeout)
        try:
            self.socket.connect((self.host, self.port))
            self.connected = True
        except socket.error as e:
            print(e, file=sys.stderr)
            if e[0] == 61:
                print("Try opening IGV and ensure port is set to %d" % self.port, file=sys.stderr)
            self.connected = False
    
    def _test(self):
        if self.send('echo', check_response='echo'):
            if not self.quiet:
                print("Test successful", file=sys.stderr)
        else:
                print("ERROR: Test failed", file=sys.stderr)
    
    def close(self):
        self.connected = False
        self.socket.close()
    
    def send(self, command, check_response='OK', attempts=5):
        if not self.connected:
            print("ERROR: Not connected to IGV.", file=sys.stderr)
            return False

        if not isinstance(command, str):
            raise TypeError(f'expected command to be str: {command}')

        success = False
        attempt = 0

        command = command.rstrip('\n') + '\n'

        while success is False:
            attempt += 1
            try:
                self.socket.send(command.encode('utf-8'))
                response = self.socket.recv(2048).decode('utf-8').rstrip('\n')
                # print(f'raw response as str: {response}', file=sys.stderr)
                if response == check_response:
                    success = True
                else:
                    print('WARNING: IGV responded "%s"' % response.strip(), file=sys.stderr)
            except socket.timeout:
                print("WARNING: socket timed out on attempt %d" % attempt, file=sys.stderr)
                self._connect_socket()
            
            if attempt >= attempts:
                break
        
        if not success:
            print("ERROR: Send failed after %d attempts." % attempt, file=sys.stderr)
        
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
        self.send("exit", confirm=False)
        self.close()
    
    def goto(self,locus):
        ''' Scrolls to a single locus or space-delimited list of loci. If a list is provided,
            these loci will be displayed in a split screen view. Use any syntax that is valid
            in the IGV search box.
        '''
        self.send("goto %s" % locus)
        return self
    
    def region(self, chr, start, end):
        ''' Not implemented '''
        pass
  
    def maxPanelHeight(self, height):
        ''' Sets the number of vertical pixels (height) of each panel to include in image
        '''
        self.send("maxPanelHeight %d" % height)
        return self
          
    def snapshot(self, filename=None):
        ''' Saves a snapshot of the IGV window to an image file. If filename is omitted,
            writes a PNG file with a filename generated based on the locus. If filename is
            specified, the filename extension determines the image file format, which must be
            .png, .jpg, or .svg
        '''
        extensions = set(['png','svg','jpg'])
        if filename is None:
            self.send('snapshot')
        else:
            if filename.split('.')[-1] not in extensions:
                raise ValueError(f'ERROR: Filename "{filename}" is invalid')
            self.send('snapshot %s' % filename)
        return self
  
    def snapshotDirectory(self,dname='.'):
        ''' Sets the directory in which to write images
        '''
        import os
        assert os.path.isdir(dname)
        self.send('snapshotDirectory %s' % os.path.abspath(dname))
        return self
    
    def viewaspairs(self, trackName=None):
        ''' Set the display mode for an alignment track to "View as pairs". trackName is
            optional.
        '''
        if trackName is None:
            self.send('viewaspairs')
        else:
            self.send('viewaspairs %s' % trackName)
        return self
  
    def squish(self, trackName=None):
        ''' Squish a given trackName. trackName is optional, and if it is not supplied all
            annotation tracks are squished.
        '''
        if trackName is None:
            self.send('squish')
        else:
            self.send('squish %s' % trackName)
        return self
  
    def collapse(self, trackName=None):
        ''' Collapse a given trackName. trackName is optional, and if it is not supplied all
            annotation tracks are collapsed.
        '''
        if trackName is None:
            self.send('collapse')
        else:
            self.send('collapse %s' % trackName)
        return self
    
    def expand(self, trackName=None):
        ''' Expand a given trackName. trackName is optional, and if it is not supplied all
            annotation tracks are expanded.
        '''
        if trackName is None:
            self.send('expand')
        else:
            self.send('expand %s' % trackName)
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
                print('Setting sleep interval to %d' % int(ms), file=sys.stderr)
        return self


def igv_init(genome, tracks=[]):
    igv = IGV()
    if not igv.connected:
        return None

    igv.new().genome(genome)
    for t in tracks:
        igv.load(t)
    return igv
