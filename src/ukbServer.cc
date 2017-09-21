#ifdef UKB_SERVER

#include "ukbServer.h"
// for deamon
#include <iostream>
#include <syslog.h>
#include <unistd.h>
#include <boost/asio/signal_set.hpp>

namespace ukb {

	int start_daemon(unsigned int port, void (*load_kb_dict)(bool), bool (*func)(sSession &)) {

		// Code "borrowed" from asio daemon example (boost license).

		int fd[2]; // Pipe for communicating btw. parent and (second) child
		int nbytes;
		char buffer[256];

		boost::asio::io_service io_service;
		sServer *server;
		try {
			// Initialise the server before becoming a daemon. If the process is
			// started from a shell, this means any errors will be reported back to the
			// user.
			server = new sServer(io_service, port, func);

		} catch (std::exception& e) {
			std::cerr << "[E] can not start daemon: " << e.what() << std::endl;
			exit(1);
		}

		if (pipe(fd)) {
			std::cerr << "[E] start_daemon: pipe failed" << std::endl;
			exit(1);
		}

		try {

			// Register signal handlers so that the daemon may be shut down. You may
			// also want to register for other signals, such as SIGHUP to trigger a
			// re-read of a configuration file.

			boost::asio::signal_set signals(io_service, SIGINT, SIGTERM, SIGHUP);
			signals.async_wait(boost::bind(&boost::asio::io_service::stop, &io_service));

			// Inform the io_service that we are about to become a daemon. The
			// io_service cleans up any internal resources, such as threads, that may
			// interfere with forking.
			io_service.notify_fork(boost::asio::io_service::fork_prepare);

			// Fork the process and have the parent exit. If the process was started
			// from a shell, this returns control to the user. Forking a new process is
			// also a prerequisite for the subsequent call to setsid().
			if (pid_t pid = fork()) {
				if (pid > 0) {
					// We're in the parent process and need to exit.
					//
					// When the exit() function is used, the program terminates without
					// invoking local variables' destructors. Only global variables are
					// destroyed. As the io_service object is a local variable, this means
					// we do not have to call:
					//
					//   io_service.notify_fork(boost::asio::io_service::fork_parent);
					//
					// However, this line should be added before each call to exit() if
					// using a global io_service object. An additional call:
					//
					//   io_service.notify_fork(boost::asio::io_service::fork_prepare);
					//
					// should also precede the second fork().
					close(fd[1]); // Parent process closes up output side of pipe
					nbytes = read(fd[0], buffer, sizeof(buffer)); // wait until receive notification from child
					if (nbytes < 0 or strcmp("OK",buffer)) {
						std::cerr << "Error initializing server\n" << buffer << "\n";
						exit(1);
					}
					exit(0);
				} else {
					syslog(LOG_ERR | LOG_USER, "First fork failed: %m");
					exit(1);
				}
			}

			close(fd[0]); // Child process closes up input side of pipe

			// Make the process a new session leader. This detaches it from the
			// terminal.
			setsid();

			// A process inherits its working directory from its parent. This could be
			// on a mounted filesystem, which means that the running daemon would
			// prevent this filesystem from being unmounted. Changing to the root
			// directory avoids this problem.
			if(chdir("/")) {
				std::cerr << "[E] start_daemon: can not chdir to /" << std::endl;
				exit(1);
			}

			// The file mode creation mask is also inherited from the parent process.
			// We don't want to restrict the permissions on files created by the
			// daemon, so the mask is cleared.
			umask(0);

			// A second fork ensures the process cannot acquire a controlling terminal.
			if (pid_t pid = fork()) {
				if (pid > 0) {
					exit(0);
				} else {
					syslog(LOG_ERR | LOG_USER, "Second fork failed: %m");
					exit(1);
				}
			}

			// Close the standard streams. This decouples the daemon from the terminal
			// that started it.
			close(0);
			close(1);
			close(2);

			// We don't want the daemon to have any standard input.
			if (open("/dev/null", O_RDONLY) < 0) {
				syslog(LOG_ERR | LOG_USER, "Unable to open /dev/null: %m");
				exit(1);
			}

			// Send standard output to a log file.
			// const char* output = "/tmp/asio.daemon.out";
			// const int flags = O_WRONLY | O_CREAT | O_APPEND;
			// const mode_t mode = S_IRUSR | S_IWUSR | S_IRGRP | S_IROTH;
			// if (open(output, flags, mode) < 0) {
			//   syslog(LOG_ERR | LOG_USER, "Unable to open output file %s: %m", output);
			//   return 1;
			// }

			// Also send standard error to the same log file.
			// if (dup(1) < 0) {
			//   syslog(LOG_ERR | LOG_USER, "Unable to dup output descriptor: %m");
			//   return 1;
			// }

			// Inform the io_service that we have finished becoming a daemon. The
			// io_service uses this opportunity to create any internal file descriptors
			// that need to be private to the new process.
			io_service.notify_fork(boost::asio::io_service::fork_child);

			// The io_service can now be used normally.

			// Pipe for communicating with parent
			// Call load callback
			try {
				(*load_kb_dict)(true); // 'true' because we call it from the daemon
			} catch (std::exception& e) {
				strncpy(buffer,e.what(), 255);
				buffer[254] = '\0'; // prevent exception messages with more than 256 chars
				if (!write(fd[1], buffer, (strlen(buffer)+1))) {
					syslog(LOG_ERR, "Error on write");
				}
				throw e;
			}
			strcpy(buffer, "OK");
			if(!write(fd[1], buffer, (strlen(buffer)+1))) {
				syslog(LOG_ERR, "Error on write. Exiting.");
				exit(1);
			}
			syslog(LOG_INFO | LOG_USER, "UKB daemon started");
			close(fd[1]); // close up also the ouput pipe
			io_service.run();
			syslog(LOG_INFO | LOG_USER, "UKB daemon stopped");
		} catch (std::exception& e) {
			syslog(LOG_ERR | LOG_USER, "[E] daemon: %s", e.what());
		}
		delete server;
		return 0;
	}

	/////////////////////////////////////////////////////////////
	// sSession class

	sSession::sSession(boost::asio::io_service& io, bool (*func)(sSession &))
		: m_func(func),
		  m_socket(io),
		  m_reader(m_socket),
		  m_writer(m_socket) {
	}

	boost::asio::ip::tcp::tcp::socket & sSession::socket() {
		return m_socket;
	}

	bool sSession::start() {
		// Start synchronous operations
		// By now, just call m_func
		return (*m_func)(*this);
	}

	bool sSession::receive(std::string & result) {
		return m_reader.read_string(result);
	}

	void sSession::send(const std::string & line) {
		m_writer.write_string(line);
	}

	//////////////////////////////////////////////////////////////
	// sServer class

	sServer::sServer(boost::asio::io_service & io, unsigned int port, bool (*func)(sSession &)) :
		m_io(io),
		m_func(func),
		m_endpoint(boost::asio::ip::tcp::tcp::v4(), port),
		m_acceptor(m_io, m_endpoint) {
		start_accept();
	}

	void sServer::start_accept() {
		sSession *new_session = new sSession(m_io, m_func);
		m_acceptor.async_accept(new_session->socket(),
								boost::bind(&sServer::handle_accept, this, new_session,
											boost::asio::placeholders::error));
	}

	void sServer::handle_accept(sSession *new_session,
								const boost::system::error_code& error) {
		bool status = true;
		if (!error) {
			status = new_session->start(); // false means finish
		}
		delete new_session;
		if (!status)
			m_io.stop();
		else
			start_accept();
	}

	//////////////////////////////////////////////////////////////
	// sClient

	sClient::sClient(const std::string & host, unsigned int port) :
		m_resolver(m_io),
		m_query(host, boost::lexical_cast<std::string>(port)), //ask the dns for this resolver
		m_socket(m_io),
		m_reader(m_socket),
		m_writer(m_socket) {

		boost::asio::ip::tcp::tcp::resolver::iterator end;

		m_error = boost::asio::error::host_not_found;
		m_endpoint_iterator = m_resolver.resolve(m_query);
		while (m_endpoint_iterator != end) {
			m_socket.close();
			m_socket.connect(*m_endpoint_iterator, m_error);
			if (!m_error) break;
			m_endpoint_iterator++;
		}
	}

	boost::system::error_code sClient::error() const { return m_error; }
	std::string sClient::error_str() const { return m_error.message(); }

	bool sClient::receive(std::string & result) {
		if (this->error())
			throw std::runtime_error(std::string("[E] sClient receive: ") + this->error_str());
		return m_reader.read_string(result);
	}

	void sClient::send(const std::string & line) {
		if (this->error())
			throw std::runtime_error(std::string("[E] sClient send: ") + this->error_str());
		m_writer.write_string(line);
	}

	//////////////////////////////////////////////////////////////
	// protocol classes

	namespace protocol {

		//////////////////////////////////////////////////////////////
		// sRead class

		sRead::sRead(boost::asio::ip::tcp::tcp::socket & socket) :
			m_socket(socket),
			m_left(&m_buf[0]),
			m_right(m_left) {}


		bool sRead::read_string(std::string & out) {
			std::string().swap(out); // empty line
			size_t N = read_size(); // get string size
			if (!N) return false; // EOF
			if (read_nstring(out, N)) {
				// Connection unexpectedly closed by peer.
				throw std::runtime_error(std::string("[E] read_string: Connection unexpectedly closed by peer.\n"));
			}
			return true;
		}

		bool sRead::read_packet() {
			size_t len = m_socket.read_some(boost::asio::buffer(m_buf), m_error); //read data
			if (m_error == boost::asio::error::eof)
				// Connection closed by peer.
				return false;
			else if (m_error)
				throw std::runtime_error(std::string("Socket error: ") + m_error.message()); // Some other error.
			m_left = &m_buf[0];
			m_right = m_left + len;
			return true;
		}

		size_t sRead::read_nstring(std::string & out, size_t N) {
			size_t k = 0;
			std::string(N, ' ').swap(out);
			while(1) {
				while(N and m_left < m_right) {
					out[k++] = *m_left++;
					--N;
				}
				if (!N) break;
				if (!read_packet()) break; // EOF
			}
			return N; // number of characters left
		}

		// A return value of zero size means EOF (since we don't allow zero-sized strings to be sent around)
		size_t sRead::read_size() {
			std::string aux;
			// first 1 byte with the size of the size string
			if (read_nstring(aux, 1) == 1) return 0; // If can't be read it is EOF
			size_t len = boost::lexical_cast<size_t>(aux);
			// read the size
			if (read_nstring(aux, len)) {
				// Connection unexpectedly closed by peer.
				throw std::runtime_error(std::string("[E] read_string: Connection unexpectedly closed by peer.\n"));
			}
			return boost::lexical_cast<size_t>(aux);
		}

		//////////////////////////////////////////////////////////////
		// sWrite class

		sWrite::sWrite(boost::asio::ip::tcp::tcp::socket & socket) :
			m_socket(socket) {}

		void sWrite::write_string(const std::string & line) {
			size_t N = line.size();
			if(!N) return;
			this->write_size(N);               // send string size
			this->write_data(line.c_str(), N); // send string
		}

		void sWrite::write_data(const char *buff, size_t len) {
			boost::system::error_code error;
			boost::asio::write(m_socket, boost::asio::buffer(buff,len),
							   boost::asio::transfer_all(), error); //send
			if (error)
				throw std::runtime_error(std::string("[E] sWrite send_data: ") + error.message());
		}

		void sWrite::write_size(size_t N) {
			std::string str = boost::lexical_cast<std::string>(N);
			size_t m = str.size();
			std::string size = boost::lexical_cast<std::string>(m);
			write_data(&size[0], 1);
			write_data(str.c_str(), m);
		}
	} // end of namespace protocol
}
#endif
