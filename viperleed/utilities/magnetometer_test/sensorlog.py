IP = "192.168.178.20"
PORT = 12348

size = 256

# import libraries
import time
import matplotlib.pyplot as plt
import pandas as pd
import socket
import numpy as np

def process_data(parts, mag_data):
    for part in parts:
        try:
            timestamp, nr, bla, x, y, z = part.split(";")
            if float(x) and float(y) and float(z):
                for n, dat in enumerate((x, y, z)):
                    mag_data[n].pop(0)
                    mag_data[n].append(float(dat))
            else:
                raise
        except:
            continue
    return


# create the zmq client and listen on port 1234
s = socket.socket()

s.bind((IP,PORT))
s.listen()
conn, addr = s.accept()
#s.connect((IP, PORT))

print("Connection established")

buffer = ""
length = 200

t = np.arange(0, length)
mag_data = [[0] * length, [0] * length, [0] * length]


# create plot
plt.ion()  # <-- work in "interactive mode"
fig, ax = plt.subplots()
fig.canvas.set_window_title("Mag Data")




while True:
    plt.cla()
    data = conn.recv(size)
    buffer += data.decode()
    # split
    parts = buffer.split("\n")
    # get back last, likely incomplete item into buffer
    buffer = parts.pop()
    process_data(parts, mag_data)

    # plot all data
    ax.plot(t, mag_data[0])
    #ax.plot(t, mag_data[1])
    #ax.plot(t, np.array(mag_data[2])**2+np.array(mag_data[1])**2+np.array(mag_data[0])**2)
    # show the plot
    plt.show()
    plt.pause(0.0001)  # <-- sets the current plot until refreshed

    time.sleep(0.00001)
