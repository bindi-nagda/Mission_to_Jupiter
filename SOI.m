clear all;
clc

D_Io = 421700;
D_Europa = 671034;
D_Ganymede = 1070412;
D_Callisto = 1882709;
D_Moon = 384.4* 10^3;
D_Amalthea = 181366;


m_Io = 8931900*10^16;
m_Europa = 4800000*10^16;
m_Ganymede = 14819000*10^16;
m_Callisto = 10759000*10^16;
m_J = 1.899*10^27;
m_Moon = 73.48*10^21;
m_Earth= 5.974*10^24;
m_Amalthea = 208*10^16;

rSOI_Io = D_Io * (m_Io/m_J)^(2/5);
rSOI_Europa = D_Europa * (m_Europa/m_J)^(2/5);
rSOI_Ganymede = D_Ganymede * (m_Ganymede/m_J)^(2/5);
rSOI_Callisto = D_Callisto * (m_Callisto/m_J)^(2/5);
rSOI_Moon = D_Moon * (m_Moon/m_Earth)^(2/5);
rSOI_Amalthea = D_Amalthea * (m_Amalthea/m_J)^(2/5);
