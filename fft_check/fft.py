import numpy as np
from scipy.fft import fft
import matplotlib.pyplot as plt
from obspy import read
####


yuk=read('../sac_files/231117_081412_PA_inc2_r2.5/*YUK8*')
#
print(yuk[0].stats,'\n')
plt.rcParams.update({'font.size': 14})

yf=fft(yuk[0].data)
xf = np.fft.fftfreq(len(yuk[0].data), 1 / 20) # 20 is sample rate

# yf.argmax()
plt.plot(xf, np.abs(yf),c='seagreen',lw=.75)
plt.title('FFT of trace for station YUK8')
plt.xlabel('Frequency (Hz)')
plt.ylabel('Amplitude')
plt.xlim(0,1)
plt.grid()
# plt.show()
plt.savefig('fft_yuk8.png',dpi=300,bbox_inches='tight', pad_inches=0.1)
