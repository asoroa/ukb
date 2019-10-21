package ukbProtocol ;

use strict;
use Encode;
use Carp qw(croak);
use Socket;     # This defines PF_INET and SOCK_STREAM

sub new {
	my $that = shift;
	my $class = ref($that) || $that;

	croak "Error: must pass port"
	  unless @_;
	my $self = {
				socket => undef,
				buff_size => 16384,
				buffer => "",
				a => 0,
				m => 0,
				reader => 0,
				writer => 0
			   };
	bless $self, $class;
	$self->_init($_[0]);
	return $self;
}

sub _read_packet {
	my $self = shift;
	my $len = sysread($self->{socket}, $self->{buffer}, $self->{buff_size});
	$self->{a} = 0;
	return undef unless defined $len;
	$self->{m} = $len;
	return 0 if not $len;
	return $len;
}

sub _read_nstring {
	my $self = shift;
	my $N = shift;
	my $str = "";
	while(1) {
		my $m = $self->{m} - $self->{a};
		if ($m) {
			$m = $N if $N < $m;
			$str .= substr($self->{buffer}, $self->{a}, $m);
			$self->{a} += $m;
			$N -= $m;
		}
		last if not $N;
		my $ok = $self->_read_packet();
		return (undef, 0) unless defined $ok;
		last if not $ok;
	}
	return ($str, $N); # (string, number of characters left)
}

sub _read_size {
	my $self = shift;
	my ($n, $aux, $len);
	($len, $n) = $self->_read_nstring(1);
	return 0 unless defined $len;
	return 0 if $n == 1;
	($aux, $n) = $self->_read_nstring($len);
	if ($n) {
		die "[E] read_size: Connection unexpectedly closed by peer.\n";
	}
	return $aux;
}

sub _read_string {
	my $self = shift;
	my $N = $self->_read_size();
	return 0 unless $N;
	my ($str, $m) = $self->_read_nstring($N);
	die "[E] Connection unexpectedly closed by peer\n" if not defined $str or $m;
	return $str;
}

sub receive {
	my $self = shift;
	return $self->_read_string();
}

sub _write_data {
	my $self = shift;
	my $buff = shift;
	my $len = shift;

	my $n = 0;
	my $off = 0;
	while(1) {
		$n = syswrite($self->{socket}, $buff, $len, $off);
		die "[E] write_data: Can't send data\n" unless $n;
		$len -= $n;
		$off += $n;
		last unless $len;
	}
}

sub _write_size {
	my $self = shift;
	my $N = shift;
	my $len = length($N);
	$self->_write_data($len, 1);
	$self->_write_data($N, $len);
}

sub send {
	my $self = shift;
	my $line = shift;
	$line = encode('UTF-8', $line) if utf8::is_utf8($line);
	my $N = length($line);
	return unless $N;
	$self->_write_size($N);
	$self->_write_data($line, $N);
}

sub _init {
	my $self = shift;
	my $port = shift;

	my $server = "localhost";  # Host IP running the server

	# create the socket, connect to the port
	socket(my $socket, PF_INET, SOCK_STREAM, (getprotobyname('tcp'))[2])
	  or die "[E] Can't create a socket $!\n";
	connect( $socket, pack_sockaddr_in($port, inet_aton($server)))
	  or die "[E] Can't connect to port $port! \n";
	binmode($socket);
	$self->{socket} = $socket;
}

(1);
