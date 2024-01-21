# About the project
<p align="justify">
The goal of the project is implementation of audio (guitar) effects which can be applied to audio signals. Project contains implementations of four effects: equalizer, wah-wah, tremolo and flanger. All of these effects are implemented using both python (python subdirectory) and C (C subdirectory). Python program can be run on any personal computer. C program is written for the targeted platform, ADSP-21489 but can be compiled with C compiler on a personal computer with some modifications. Project is made as a part of a grade for university subject "Systems for digital signal processing".
</p>

## How to run the programs
### Python
<p align="justify">
Python implementation is contained in a jupyter notebook file "Efekti.ipynb". To run this file you need to install python kernel with NumPy, SciPy and Matplotlib libraries. After that you can load the file into enviroment like Jupyter, Spider or VS Code and run the cells.
</p>

### C
<p align="justify">
To run the C program you need to install CrossCore Embedded Studio (CCES) development enviroment and open the C subfolder as a project inside it. Apart from the CCES you need ADSP-21489 EzKit board as algorithm selection and parameter configuration is done by pushing buttons available on the platform. After loading the project and connecting the board run the application by clicking Debug.
</p>

## Want to test your own audio signals ?
<p align="justify">
You can load your signals into the notebook file by using any python integrated function to read audio files (for example .wav files you can read using scipy wavfile.read function). Be careful to load sample rate and audio samples to sample_freq variable and compound_signal array respectively (Check the examples where acoustic.wav is being loaded).
To apply effects in CCES project you need to export raw audio signal samples into C array and define it's length into a C header file compound_signal.h. Be aware of memory limitations...
</p>
