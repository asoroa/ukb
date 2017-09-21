// -*-C++-*-

#ifndef UKB_SERVER_HPP
#define UKB_SERVER_HPP

#ifdef UKB_SERVER

#include <string>
#include <boost/asio.hpp>
#include <boost/bind.hpp>
#include <boost/lexical_cast.hpp>

// Class to handle client/server operation in ukb

namespace ukb {

	class sSession;

	int start_daemon(unsigned int port, void (*pre)(bool), bool (*func)(sSession &));

	namespace protocol {
		// Protocol is as follows:
		//
		// - one byte (M).
		// - a string of length M which describes an unsigned int (N).
		// - a string of length N with the actual data.
		//
		// for instance, to send/receive the string "This is a string" (whose lenght
		// is 16)
		//
		// - first, a byte with value 2
		// - then a string with value "16"
		// - finally, a string with value "This is a string"

		// sRead
		// read strings from socket following protocol
		class sRead {

		public:
			sRead(boost::asio::ip::tcp::tcp::socket & socket);
			bool read_string(std::string & out);

		private:
			bool read_packet();
			size_t read_nstring(std::string & out, size_t N);
			// A return value of zero size means EOF (since we don't allow zero-sized strings to be sent around)
			size_t read_size();

			static const unsigned int buff_size  = 16384; //size of the send buffer
			char m_buf[buff_size];
			boost::asio::ip::tcp::tcp::socket & m_socket;
			boost::system::error_code m_error;
			char *m_left, *m_right; // actual range

		};

		// sWrite
		// write strings from socket following protocol
		class sWrite {

		public:

			sWrite(boost::asio::ip::tcp::tcp::socket & socket);
			void write_string(const std::string & line);

		private:
			boost::asio::ip::tcp::tcp::socket & m_socket;

			void write_data(const char *buff, size_t len);
			void write_size(size_t N);
		};
	}

	// A server session. Create socket and communicate through it.
	//
	// func is a callback function responsible of dealing with server side
	// communication. It is typically defined inside the main file (ukb_wsd or
	// ukb_ppv)

	class sSession {

	public:
		sSession(boost::asio::io_service& io, bool (*func)(sSession &));

		boost::asio::ip::tcp::tcp::socket & socket();
		bool start();
		bool receive(std::string & result);
		void send(const std::string & line);

	private:

		bool (*m_func)(sSession &);
		boost::asio::ip::tcp::tcp::socket m_socket;
		protocol::sRead m_reader;
		protocol::sWrite m_writer;

	};

	// Server main class. Accept connections asyncronously, create a session and
	// launch it.

	class sServer {

	public:
		sServer(boost::asio::io_service & io, unsigned int port, bool (*func)(sSession &));

	private:

		void start_accept();
		void handle_accept(sSession *new_session,
						   const boost::system::error_code& error);

		// Connection stuff
		boost::asio::io_service & m_io; //main asio object
		bool (*m_func)(sSession &);
		boost::asio::ip::tcp::tcp::endpoint m_endpoint;
		boost::asio::ip::tcp::tcp::acceptor m_acceptor;

	};

	// Client side of communication.
	//
	// create a session with the server (blocking until communication is done) and
	// offer functions for sendig/receiving strings.

	class sClient {

	public:
		// Connects to host:port. Sets error on failure.
		sClient(const std::string & host, unsigned int port);

		boost::system::error_code error() const;
		std::string error_str() const;

		bool receive(std::string & result);
		void send(const std::string & line);
	private:

		boost::system::error_code m_error; // zero if connection is succesful
		boost::asio::io_service m_io; //asio main object
		boost::asio::ip::tcp::tcp::resolver m_resolver;
		boost::asio::ip::tcp::tcp::resolver::query m_query;
		boost::asio::ip::tcp::tcp::resolver::iterator m_endpoint_iterator; //iterator if multiple answers for a given name
		boost::asio::ip::tcp::tcp::socket m_socket;

		protocol::sRead m_reader;
		protocol::sWrite m_writer;
	};
}

#endif // UKB_SERVER

#endif // UKB_SERVER_HPP
