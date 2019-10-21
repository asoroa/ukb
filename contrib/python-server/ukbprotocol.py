#!/usr/bin/env python3

import socket

class UkbSession:
    def __init__(self, port, server = "localhost"):
        self.buffer_size = 16384
        self.buffer = None
        self.a = 0
        self.b = 0
        try:
            self.socket = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
            self.socket.connect((socket.gethostbyname(server), port))
        except:
            raise OSError("Can not connect to server '{}:{}'".format(server, str(port)))

    def _read_packet(self):
        '''Reads a packet up to self.buffer_size and leaves it in self.buffer. Returns the number of bytes read.'''
        self.buffer = self.socket.recv(self.buffer_size)
        if not self.buffer:
            raise OSError
        self.a = 0
        self.b = len(self.buffer)

    def _read_nstring(self, N):
        '''Tries to read up to N bytes from buffer. Return the number of bytes left.'''
        mydata = b''
        while True:
            m = min(self.b - self.a, N)
            if m:
                mydata += self.buffer[self.a:self.a + m]
                self.a += m
                N -= m
            if N == 0:
                break
            self._read_packet()
        return mydata

    def _read_size(self):
        '''Read the size of the incoming string. Returns an int'''
        len_bytes = self._read_nstring(1)
        size = 0
        try:
            aux = self._read_nstring(int(len_bytes))
            size = int(aux)
        except ValueError:
            raise Exception('_read_size: Protocol error')
        return size

    def _write_data(self, data):
        '''Send data bytes'''
        N = len(data)
        a = 0
        while a < N:
            n = self.socket.send(data[a:])
            if not n:
                raise OSError
            a += n

    def _write_size(self, N):
        '''Send an int following protocol.'''
        length_string = str(N)
        self._write_data(str(len(length_string)).encode('ascii')) # first, one byte string with the length of the length-string
        self._write_data(length_string.encode('ascii')) # then, the length-string as bytes

    def recv(self):
        '''Receive a string from the socket. Returns a utf-8 encoded string'''
        if not self.socket:
            raise OSError('recv: no socket')
        data = b''
        try:
            l = self._read_size()
            if l == 0: return ''
            data = self._read_nstring(l)
        except OSError:
            raise OSError('recv: connection closed by peer')
        return data.decode('utf-8', errors='replace')

    def send(self, data):
        '''Send a string to the server, following ukb protocol'''
        if not self.socket:
            raise OSError('send: no socket')
        if isinstance(data, str):
            data = data.encode('utf-8')
        l = len(data)
        if l == 0: return
        try:
            self._write_size(l)
            self._write_data(data)
        except OSError:
            raise OSError('send: connection closed by peer')

    def close(self):
        if self.socket:
            self.socket.close()
        self.socket = None
