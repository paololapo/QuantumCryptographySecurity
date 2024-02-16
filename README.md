# Quantum Cryptography & Security
Analysis and results of laboratory experiments for the 'Quantum Cryptography and Security' course (AY 2023/24). <br>
This repository is organized into three folders, one for each experiment: *Quantum Random Numbers Generation*, *Error Correction*, and *Quantum Key Distribution*. Each folder has a <code> functions.py </code> file with the analysis, a folder with the graphs, and the final report. 

## Introductions
### Quantum Random Numbers Generation
The implementation of QRNG can be characterized by the degree of trust in the different elements of the protocol. The simplest case is the *trusted setup*, in which all the elements are supposed to be controlled and uncorrelated with the environment. We can relax this hypothesis and consider *semi-Device-Independent* setups with uncharacterized sources or measurements. <br> 
In this experiment, we use the phenomenon known as *spontaneous parametric downconversion* to generate a two-photon entangled state and characterize it in terms of polarization. The security analysis is based on the *leftover hash lemma* with bounds on the min-entropy based on either the state tomography or the entropic uncertainty principle. 

### Error correction
To reconcile the keys, Alice and Bob need to cooperate over a classic public channel. In this report, we describe three different protocols: *Cascade* (based on parity comparison and iterative binary search), *Winnow* (based on parity comparison and syndrome decoding), and the *LDPC codes* (with rate modulation via puncturing and shortening). <br>
We compare the performances in terms of two efficiency metrics based on the error correction capability and the *Slepian-Wolf bound*.

### Quantum Key Distribution
The standard QKD protocols are designed to work with true single-photons. However, such experimental setups are still unavailable for practical implementation, and weak coherent laser pulses (vulnerable to the so-called *photon number splitting* attack) are used instead. <br>
In this report, we will introduce and analyze a *3-states 1-decoy* QKD protocol, evaluating the security of the obtained keys in the finite scenario. Such methods provide robust protocols to overcome the limitations of the multi-photon events and the finite keys.

## References
### QRNG
[1]. Giuseppe Vallone, Davide G Marangon, Marco Tomasin, Paolo Villoresi, "Quantum randomness certified by the uncertainty principle," in *Physical Review A*, vol. 90, no. 5, pp. 052327, APS, 2014.

[2]. Joseph B Altepeter, Evan R Jeffrey, Paul G Kwiat, "Photonic state tomography," in *Advances in Atomic, Molecular, and Optical Physics*, vol. 52, pp. 105–159, Elsevier, 2005.

[3]. Daniel FV James, Paul G Kwiat, William J Munro, Andrew G White, "Measurement of qubits," in *Physical Review A*, vol. 64, no. 5, pp. 052312, APS, 2001.

[4]. M Fiorentino, C Santori, SM Spillane, RG Beausoleil, WJ Munro, "Secure self-calibrating quantum random-bit generator," in *Physical Review A*, vol. 75, no. 3, pp. 032334, APS, 2007.

[5]. Xiongfeng Ma, Feihu Xu, He Xu, Xiaoqing Tan, Bing Qi, Hoi-Kwong Lo, "Postprocessing for quantum random-number generators: Entropy evaluation and randomness extraction," in *Physical Review A*, vol. 87, no. 6, pp. 062327, APS, 2013.

[6]. Marco Tomamichel, Christian Schaffner, Adam Smith, Renato Renner, "Leftover hashing against quantum side information," in *IEEE Transactions on Information Theory*, vol. 57, no. 8, pp. 5524–5535, IEEE, 2011.


### EC
[1]. Miralem Mehic, Marcin Niemiec, Harun Siljak, Miroslav Voznak, "Error Reconciliation in Quantum Key Distribution Protocols," 2020.

[2]. Gilles Brassard and Louis Salvail, "Secret-key reconciliation by public discussion," in *Workshop on the Theory and Application of Cryptographic Techniques*, pp. 410–423, Springer, 1993.

[3]. Tomohiro Sugimoto and Kouichi Yamazaki, "A study on secret key reconciliation protocol," in *IEICE Transactions on Fundamentals of Electronics, Communications and Computer Sciences*, vol. 83, no. 10, pp. 1987–1991, The Institute of Electronics, Information and Communication Engineers, 2000.

[4]. William T Buttler, Steven K Lamoreaux, Justin R Torgerson, GH Nickel, CH Donahue, Charles G Peterson, "Fast, efficient error reconciliation for quantum cryptography," in *Physical Review A*, vol. 67, no. 5, pp. 052303, APS, 2003.

[5]. David Elkouss, Jesus Martinez-Mateo, Vicente Martin, "Information reconciliation for quantum key distribution," *arXiv preprint arXiv:1007.1616*, 2010.

[6]. David Slepian and Jack Wolf, "Noiseless coding of correlated information sources," in *IEEE Transactions on Information Theory*, vol. 19, no. 4, pp. 471–480, IEEE, 1973.

### QKD
[1]. C. H. Bennett and G. Brassard, "Quantum Cryptography: Public Key Distribution and Coin Tossing," in *Proceedings of IEEE International Conference on Computers, Systems and Signal Processing*, vol. 175, pp. 8, New York, 1984.

[2]. AA Gaidash, VI Egorov, AV Gleim, "Revealing of photon-number splitting attack on quantum key distribution system by photon-number resolving devices," in *Journal of Physics: Conference Series*, vol. 735, no. 1, pp. 012072, IOP Publishing, 2016.

[3]. Davide Rusca, Alberto Boaron, Fadri Gr{\"u}nenfelder, Anthony Martin, Hugo Zbinden, "Finite-key analysis for the 1-decoy state QKD protocol," in *Applied Physics Letters*, vol. 112, no. 17, AIP Publishing, 2018.

[4]. Charles Ci Wen Lim, Marcos Curty, Nino Walenta, Feihu Xu, Hugo Zbinden, "Concise security bounds for practical decoy-state quantum key distribution," in *Physical Review A*, vol. 89, no. 2, pp. 022307, APS, 2014.

[5]. Wassily Hoeffding, "Probability Inequalities for Sums of Bounded Random Variables," in *Journal of the American Statistical Association*, vol. 58, no. 301, pp. 13–30, [American Statistical Association, Taylor & Francis, Ltd.], 1963.

[6]. Costantino Agnesi, Marco Avesani, Andrea Stanco, Paolo Villoresi, Giuseppe Vallone, "All-fiber self-compensating polarization encoder for quantum key distribution," in *Optics letters*, vol. 44, no. 10, pp. 2398–2401, Optica Publishing Group, 2019.

[7]. Marco Avesani, Luca Calderaro, Giulio Foletto, Costantino Agnesi, Francesco Picciariello, Francesco BL Santagiustina, Alessia Scriminich, Andrea Stanco, Francesco Vedovato, Mujtaba Zahidy, "Resource-effective quantum key distribution: a field trial in Padua city center," in *Optics letters*, vol. 46, no. 12, pp. 2848–2851, Optica Publishing Group, 2021.
