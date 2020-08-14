import http.server
import socketserver
from threading import Timer
import webbrowser
"""
This script finds the first open port between 8000 and 8100
and starts a http server and then opens a new browser tab with localhost
Code from:
https://docs.python.org/3/library/webbrowser.html
https://codereview.stackexchange.com/questions/116450/find-available-ports-on-localhost
https://docs.python.org/3/library/http.server.html

"""

def is_port_in_use(port):
    import socket
    with socket.socket(socket.AF_INET, socket.SOCK_STREAM) as s:
        return s.connect_ex(('localhost', port)) == 0

def check_port():
    for i in range(8000, 8100):
            if is_port_in_use(i) == False:
                port = i
                break
    return i

def run_server(PORT):
    Handler = http.server.SimpleHTTPRequestHandler
    with socketserver.TCPServer(("", PORT), Handler) as httpd:
        print("serving at port", PORT)
        httpd.serve_forever()

def open_localhost(port, path):
    webbrowser.open_new(f"http://localhost:{port}/{path}")

def run_sunburst(path='sunburst/'):
    port = check_port()
    Timer(1, open_localhost, args=[port, path]).start();
    run_server(port)

if __name__ == "__main__":
    run_sunburst()
