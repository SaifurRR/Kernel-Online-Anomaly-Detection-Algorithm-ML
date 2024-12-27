# Kernel-Online-Anomaly-Detection-Algorithm
Contemporary patient surveillance systems have streamlined central surveillance into the electronic health record interface. They are able to process the sheer volume of patient data by adopting machine learning approaches. However, these systems are not suitable for implementation in many hospitals, mostly in developing countries, with limited human, financial, and technological resources. Through conducting thorough research on intensive care facilities, we designed a novel central patient monitoring system and in this paper, we describe the working prototype of our system. The proposed prototype comprises inexpensive peripherals and a simplistic user interface. Our central patient monitoring system implements Kernel-based On-line Anomaly Detection (KOAD) algorithm for emergency event signaling. By evaluating continuous patient data, we show that the system is able to detect critical events in real-time reliably and has a low false alarm rate.


#### 1. Central Station of Automated Signaling System in Intensive Care Unit (ICU)
  
![ALT TEXT](https://github.com/SaifurRR/Kernel-Online-Anomaly-Detection-Algorithm-ML/blob/main/import/1_Automated_Signaling_System_ICU_Central.png)****     


#### 2. Waveform Collection List from Philips IntelliVue MX750 
![ALT TEXT](https://github.com/SaifurRR/Kernel-Online-Anomaly-Detection-Algorithm-ML/blob/main/import/2_Waveforms_collected_from_Patient_Monitoring_Unit.png)****


#### 3. Data Export Test Tool (DETT) for Philips IntelliVue Patient Monitoring Unit (PMU)
 
![ALT TEXT](https://github.com/SaifurRR/Kernel-Online-Anomaly-Detection-Algorithm-ML/blob/main/import/3_Data_Extraction_Tool_Philips.png)****

#### 4. Data Export Test Tool (DETT) data format from Philips IntelliVue Patient Monitoring Unit (PMU)
 
![ALT TEXT](https://github.com/SaifurRR/Kernel-Online-Anomaly-Detection-Algorithm-ML/blob/main/import/5_DETT_data_%20format.png)****

#### 5. Kernel Online Anomaly Detection (KOAD) Algorithm 
 
![ALT TEXT](https://github.com/SaifurRR/Kernel-Online-Anomaly-Detection-Algorithm-ML/blob/main/import/4_Detection_Statistics_Koad_Algorithm.png)****

Progression of the detection statistic δt, with KOAD run over a set of 8 ion concentration levels of a patient across 100 timesteps. Timesteps corresponding to the identified critical events are indicated as red, filled stems. Sample experiment with ν1 = 0.05, ν2 = 0.20
