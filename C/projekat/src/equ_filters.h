/** @file equ_filters.h
 *
 * 	@brief Header file containing 5 band equalizer FIR filters
 */

#define NUM_TAPS 500

/*
 * Low pass filter coefficients with cutoff frequency of 300Hz
 */
const float pm FILTER_1[NUM_TAPS] = {
	-0.0008461,-0.0002454,-0.0002719,-0.0002930,-0.0003077,
	-0.0003144,-0.0003129,-0.0003021,-0.0002820,-0.0002516,
	-0.0002116,-0.0001625,-0.0001072,-0.0000442,0.0000239,
	0.0000934,0.0001647,0.0002336,0.0002991,0.0003577,
	0.0004079,0.0004472,0.0004746,0.0004882,0.0004873,
	0.0004716,0.0004420,0.0003980,0.0003420,0.0002753,
	0.0002004,0.0001200,0.0000371,-0.0000451,-0.0001233,
	-0.0001943,-0.0002551,-0.0003030,-0.0003354,-0.0003510,
	-0.0003484,-0.0003272,-0.0002879,-0.0002314,-0.0001598,
	-0.0000754,0.0000185,0.0001185,0.0002203,0.0003199,
	0.0004130,0.0004955,0.0005634,0.0006136,0.0006431,
	0.0006500,0.0006333,0.0005925,0.0005287,0.0004436,
	0.0003399,0.0002211,0.0000917,-0.0000435,-0.0001791,
	-0.0003096,-0.0004296,-0.0005337,-0.0006172,-0.0006761,
	-0.0007071,-0.0007081,-0.0006782,-0.0006175,-0.0005277,
	-0.0004116,-0.0002730,-0.0001171,0.0000502,0.0002225,
	0.0003927,0.0005537,0.0006985,0.0008208,0.0009145,
	0.0009750,0.0009987,0.0009834,0.0009285,0.0008350,
	0.0007055,0.0005442,0.0003567,0.0001499,-0.0000683,
	-0.0002893,-0.0005040,-0.0007036,-0.0008793,-0.0010232,
	-0.0011284,-0.0011895,-0.0012025,-0.0011655,-0.0010786,
	-0.0009439,-0.0007654,-0.0005492,-0.0003031,-0.0000364,
	0.0002405,0.0005166,0.0007805,0.0010208,0.0012270,
	0.0013895,0.0015003,0.0015533,0.0015446,0.0014727,
	0.0013387,0.0011461,0.0009013,0.0006126,0.0002907,
	-0.0000523,-0.0004027,-0.0007464,-0.0010691,-0.0013568,
	-0.0015967,-0.0017776,-0.0018901,-0.0019276,-0.0018865,
	-0.0017659,-0.0015686,-0.0013003,-0.0009701,-0.0005895,
	-0.0001729,0.0002637,0.0007030,0.0011269,0.0015175,
	0.0018577,0.0021320,0.0023270,0.0024324,0.0024410,
	0.0023495,0.0021585,0.0018729,0.0015015,0.0010570,
	0.0005552,0.0000151,-0.0005423,-0.0010947,-0.0016191,
	-0.0020930,-0.0024955,-0.0028075,-0.0030135,-0.0031014,
	-0.0030639,-0.0028984,-0.0026075,-0.0021990,-0.0016858,
	-0.0010857,-0.0004202,0.0002854,0.0010034,0.0017048,
	0.0023602,0.0029409,0.0034206,0.0037759,0.0039877,
	0.0040422,0.0039316,0.0036542,0.0032154,0.0026271,
	0.0019079,0.0010823,0.0001803,-0.0007643,-0.0017142,
	-0.0026307,-0.0034749,-0.0042092,-0.0047990,-0.0052141,
	-0.0054301,-0.0054297,-0.0052035,-0.0047510,-0.0040804,
	-0.0032094,-0.0021644,-0.0009800,0.0003021,0.0016341,
	0.0029641,0.0042379,0.0054005,0.0063989,0.0071837,
	0.0077113,0.0079461,0.0078619,0.0074436,0.0066883,
	0.0056061,0.0042202,0.0025669,0.0006953,-0.0013343,
	-0.0034518,-0.0055787,-0.0076309,-0.0095213,-0.0111620,
	-0.0124677,-0.0133583,-0.0137618,-0.0136172,-0.0128763,
	-0.0115064,-0.0094914,-0.0068334,-0.0035529,0.0003109,
	0.0047008,0.0095427,0.0147470,0.0202112,0.0258224,
	0.0314598,0.0369988,0.0423137,0.0472816,0.0517858,
	0.0557192,0.0589872,0.0615107,0.0632281,0.0640975,
	0.0640975,0.0632281,0.0615107,0.0589872,0.0557192,
	0.0517858,0.0472816,0.0423137,0.0369988,0.0314598,
	0.0258224,0.0202112,0.0147470,0.0095427,0.0047008,
	0.0003109,-0.0035529,-0.0068334,-0.0094914,-0.0115064,
	-0.0128763,-0.0136172,-0.0137618,-0.0133583,-0.0124677,
	-0.0111620,-0.0095213,-0.0076309,-0.0055787,-0.0034518,
	-0.0013343,0.0006953,0.0025669,0.0042202,0.0056061,
	0.0066883,0.0074436,0.0078619,0.0079461,0.0077113,
	0.0071837,0.0063989,0.0054005,0.0042379,0.0029641,
	0.0016341,0.0003021,-0.0009800,-0.0021644,-0.0032094,
	-0.0040804,-0.0047510,-0.0052035,-0.0054297,-0.0054301,
	-0.0052141,-0.0047990,-0.0042092,-0.0034749,-0.0026307,
	-0.0017142,-0.0007643,0.0001803,0.0010823,0.0019079,
	0.0026271,0.0032154,0.0036542,0.0039316,0.0040422,
	0.0039877,0.0037759,0.0034206,0.0029409,0.0023602,
	0.0017048,0.0010034,0.0002854,-0.0004202,-0.0010857,
	-0.0016858,-0.0021990,-0.0026075,-0.0028984,-0.0030639,
	-0.0031014,-0.0030135,-0.0028075,-0.0024955,-0.0020930,
	-0.0016191,-0.0010947,-0.0005423,0.0000151,0.0005552,
	0.0010570,0.0015015,0.0018729,0.0021585,0.0023495,
	0.0024410,0.0024324,0.0023270,0.0021320,0.0018577,
	0.0015175,0.0011269,0.0007030,0.0002637,-0.0001729,
	-0.0005895,-0.0009701,-0.0013003,-0.0015686,-0.0017659,
	-0.0018865,-0.0019276,-0.0018901,-0.0017776,-0.0015967,
	-0.0013568,-0.0010691,-0.0007464,-0.0004027,-0.0000523,
	0.0002907,0.0006126,0.0009013,0.0011461,0.0013387,
	0.0014727,0.0015446,0.0015533,0.0015003,0.0013895,
	0.0012270,0.0010208,0.0007805,0.0005166,0.0002405,
	-0.0000364,-0.0003031,-0.0005492,-0.0007654,-0.0009439,
	-0.0010786,-0.0011655,-0.0012025,-0.0011895,-0.0011284,
	-0.0010232,-0.0008793,-0.0007036,-0.0005040,-0.0002893,
	-0.0000683,0.0001499,0.0003567,0.0005442,0.0007055,
	0.0008350,0.0009285,0.0009834,0.0009987,0.0009750,
	0.0009145,0.0008208,0.0006985,0.0005537,0.0003927,
	0.0002225,0.0000502,-0.0001171,-0.0002730,-0.0004116,
	-0.0005277,-0.0006175,-0.0006782,-0.0007081,-0.0007071,
	-0.0006761,-0.0006172,-0.0005337,-0.0004296,-0.0003096,
	-0.0001791,-0.0000435,0.0000917,0.0002211,0.0003399,
	0.0004436,0.0005287,0.0005925,0.0006333,0.0006500,
	0.0006431,0.0006136,0.0005634,0.0004955,0.0004130,
	0.0003199,0.0002203,0.0001185,0.0000185,-0.0000754,
	-0.0001598,-0.0002314,-0.0002879,-0.0003272,-0.0003484,
	-0.0003510,-0.0003354,-0.0003030,-0.0002551,-0.0001943,
	-0.0001233,-0.0000451,0.0000371,0.0001200,0.0002004,
	0.0002753,0.0003420,0.0003980,0.0004420,0.0004716,
	0.0004873,0.0004882,0.0004746,0.0004472,0.0004079,
	0.0003577,0.0002991,0.0002336,0.0001647,0.0000934,
	0.0000239,-0.0000442,-0.0001072,-0.0001625,-0.0002116,
	-0.0002516,-0.0002820,-0.0003021,-0.0003129,-0.0003144,
	-0.0003077,-0.0002930,-0.0002719,-0.0002454,-0.0008461
};

/*
 * Band pass filter coefficients with (300 Hz - 600 Hz) bandwidth
 */
const float pm FILTER_2[NUM_TAPS] = {
	-0.0006680,0.0004082,0.0002860,0.0001884,0.0001103,
	0.0000490,0.0000031,-0.0000260,-0.0000366,-0.0000260,
	0.0000068,0.0000620,0.0001358,0.0002221,0.0003109,
	0.0003907,0.0004478,0.0004706,0.0004481,0.0003761,
	0.0002527,0.0000862,-0.0001117,-0.0003245,-0.0005270,
	-0.0007010,-0.0008243,-0.0008806,-0.0008606,-0.0007652,
	-0.0006018,-0.0003869,-0.0001426,0.0001042,0.0003269,
	0.0005021,0.0006138,0.0006537,0.0006241,0.0005357,
	0.0004082,0.0002651,0.0001332,0.0000347,-0.0000124,
	-0.0000006,0.0000676,0.0001799,0.0003146,0.0004447,
	0.0005424,0.0005821,0.0005454,0.0004244,0.0002237,
	-0.0000402,-0.0003402,-0.0006424,-0.0009099,-0.0011089,
	-0.0012133,-0.0012087,-0.0010949,-0.0008857,-0.0006079,
	-0.0002971,0.0000067,0.0002663,0.0004512,0.0005435,
	0.0005402,0.0004535,0.0003092,0.0001434,-0.0000043,
	-0.0000965,-0.0001043,-0.0000128,0.0001760,0.0004425,
	0.0007520,0.0010587,0.0013126,0.0014671,0.0014855,
	0.0013487,0.0010571,0.0006332,0.0001182,-0.0004325,
	-0.0009569,-0.0013949,-0.0016974,-0.0018326,-0.0017914,
	-0.0015886,-0.0012611,-0.0008625,-0.0004549,-0.0000997,
	0.0001521,0.0002678,0.0002387,0.0000813,-0.0001639,
	-0.0004394,-0.0006788,-0.0008178,-0.0008047,-0.0006095,
	-0.0002302,0.0003055,0.0009424,0.0016049,0.0022068,
	0.0026639,0.0029065,0.0028897,0.0026015,0.0020653,
	0.0013376,0.0005011,-0.0003469,-0.0011081,-0.0016969,
	-0.0020538,-0.0021530,-0.0020072,-0.0016655,-0.0012055,
	-0.0007220,-0.0003110,-0.0000556,-0.0000111,-0.0001956,
	-0.0005859,-0.0011193,-0.0017025,-0.0022244,-0.0025740,
	-0.0026571,-0.0024129,-0.0018249,-0.0009262,0.0002028,
	0.0014435,0.0026547,0.0036920,0.0044283,0.0047723,
	0.0046826,0.0041743,0.0033181,0.0022304,0.0010572,
	-0.0000473,-0.0009438,-0.0015282,-0.0017484,-0.0016122,
	-0.0011868,-0.0005887,0.0000345,0.0005269,0.0007501,
	0.0006066,0.0000582,-0.0008653,-0.0020660,-0.0033898,
	-0.0046470,-0.0056391,-0.0061884,-0.0061665,-0.0055151,
	-0.0042587,-0.0025050,-0.0004329,0.0017305,0.0037402,
	0.0053678,0.0064344,0.0068366,0.0065631,0.0056966,
	0.0044026,0.0029056,0.0014545,0.0002852,-0.0004176,
	-0.0005502,-0.0001098,0.0008023,0.0019926,0.0032028,
	0.0041483,0.0045631,0.0042438,0.0030864,0.0011104,
	-0.0015356,-0.0045865,-0.0076920,-0.0104625,-0.0125225,
	-0.0135661,-0.0134038,-0.0119951,-0.0094611,-0.0060740,
	-0.0022239,0.0016325,0.0050360,0.0075908,0.0090249,
	0.0092339,0.0083021,0.0064959,0.0042299,0.0020076,
	0.0003462,-0.0003061,0.0003524,0.0024154,0.0057347,
	0.0099183,0.0143658,0.0183358,0.0210406,0.0217550,
	0.0199260,0.0152694,0.0078366,-0.0019554,-0.0133413,
	-0.0252695,-0.0365053,-0.0457652,-0.0518643,-0.0538630,
	-0.0511919,-0.0437396,-0.0318916,-0.0165133,0.0011228,
	0.0194590,0.0368100,0.0515393,0.0622336,0.0678551,
	0.0678551,0.0622336,0.0515393,0.0368100,0.0194590,
	0.0011228,-0.0165133,-0.0318916,-0.0437396,-0.0511919,
	-0.0538630,-0.0518643,-0.0457652,-0.0365053,-0.0252695,
	-0.0133413,-0.0019554,0.0078366,0.0152694,0.0199260,
	0.0217550,0.0210406,0.0183358,0.0143658,0.0099183,
	0.0057347,0.0024154,0.0003524,-0.0003061,0.0003462,
	0.0020076,0.0042299,0.0064959,0.0083021,0.0092339,
	0.0090249,0.0075908,0.0050360,0.0016325,-0.0022239,
	-0.0060740,-0.0094611,-0.0119951,-0.0134038,-0.0135661,
	-0.0125225,-0.0104625,-0.0076920,-0.0045865,-0.0015356,
	0.0011104,0.0030864,0.0042438,0.0045631,0.0041483,
	0.0032028,0.0019926,0.0008023,-0.0001098,-0.0005502,
	-0.0004176,0.0002852,0.0014545,0.0029056,0.0044026,
	0.0056966,0.0065631,0.0068366,0.0064344,0.0053678,
	0.0037402,0.0017305,-0.0004329,-0.0025050,-0.0042587,
	-0.0055151,-0.0061665,-0.0061884,-0.0056391,-0.0046470,
	-0.0033898,-0.0020660,-0.0008653,0.0000582,0.0006066,
	0.0007501,0.0005269,0.0000345,-0.0005887,-0.0011868,
	-0.0016122,-0.0017484,-0.0015282,-0.0009438,-0.0000473,
	0.0010572,0.0022304,0.0033181,0.0041743,0.0046826,
	0.0047723,0.0044283,0.0036920,0.0026547,0.0014435,
	0.0002028,-0.0009262,-0.0018249,-0.0024129,-0.0026571,
	-0.0025740,-0.0022244,-0.0017025,-0.0011193,-0.0005859,
	-0.0001956,-0.0000111,-0.0000556,-0.0003110,-0.0007220,
	-0.0012055,-0.0016655,-0.0020072,-0.0021530,-0.0020538,
	-0.0016969,-0.0011081,-0.0003469,0.0005011,0.0013376,
	0.0020653,0.0026015,0.0028897,0.0029065,0.0026639,
	0.0022068,0.0016049,0.0009424,0.0003055,-0.0002302,
	-0.0006095,-0.0008047,-0.0008178,-0.0006788,-0.0004394,
	-0.0001639,0.0000813,0.0002387,0.0002678,0.0001521,
	-0.0000997,-0.0004549,-0.0008625,-0.0012611,-0.0015886,
	-0.0017914,-0.0018326,-0.0016974,-0.0013949,-0.0009569,
	-0.0004325,0.0001182,0.0006332,0.0010571,0.0013487,
	0.0014855,0.0014671,0.0013126,0.0010587,0.0007520,
	0.0004425,0.0001760,-0.0000128,-0.0001043,-0.0000965,
	-0.0000043,0.0001434,0.0003092,0.0004535,0.0005402,
	0.0005435,0.0004512,0.0002663,0.0000067,-0.0002971,
	-0.0006079,-0.0008857,-0.0010949,-0.0012087,-0.0012133,
	-0.0011089,-0.0009099,-0.0006424,-0.0003402,-0.0000402,
	0.0002237,0.0004244,0.0005454,0.0005821,0.0005424,
	0.0004447,0.0003146,0.0001799,0.0000676,-0.0000006,
	-0.0000124,0.0000347,0.0001332,0.0002651,0.0004082,
	0.0005357,0.0006241,0.0006537,0.0006138,0.0005021,
	0.0003269,0.0001042,-0.0001426,-0.0003869,-0.0006018,
	-0.0007652,-0.0008606,-0.0008806,-0.0008243,-0.0007010,
	-0.0005270,-0.0003245,-0.0001117,0.0000862,0.0002527,
	0.0003761,0.0004481,0.0004706,0.0004478,0.0003907,
	0.0003109,0.0002221,0.0001358,0.0000620,0.0000068,
	-0.0000260,-0.0000366,-0.0000260,0.0000031,0.0000490,
	0.0001103,0.0001884,0.0002860,0.0004082,-0.0006680	
};

/*
 * Band pass filter coefficients with (600 Hz - 900 Hz) bandwidth
 */
const float pm FILTER_3[NUM_TAPS] = {
	0.0004589,-0.0006407,-0.0003279,-0.0001608,-0.0000714,
	-0.0000293,-0.0000221,-0.0000436,-0.0000843,-0.0001285,
	-0.0001546,-0.0001415,-0.0000752,0.0000423,0.0001908,
	0.0003329,0.0004238,0.0004232,0.0003104,0.0000927,
	-0.0001896,-0.0004722,-0.0006802,-0.0007505,-0.0006480,
	-0.0003830,-0.0000064,0.0003947,0.0007234,0.0008985,
	0.0008704,0.0006481,0.0002852,-0.0001279,-0.0004890,
	-0.0007129,-0.0007570,-0.0006285,-0.0003807,-0.0000946,
	0.0001455,0.0002802,0.0002915,0.0002072,0.0000868,
	-0.0000004,-0.0000028,0.0000946,0.0002612,0.0004305,
	0.0005211,0.0004675,0.0002447,-0.0001153,-0.0005302,
	-0.0008855,-0.0010713,-0.0010154,-0.0007085,-0.0002128,
	0.0003541,0.0008504,0.0011506,0.0011822,0.0009471,
	0.0005196,0.0000229,-0.0004081,-0.0006690,-0.0007163,
	-0.0005774,-0.0003370,-0.0001057,0.0000214,-0.0000005,
	-0.0001489,-0.0003412,-0.0004646,-0.0004187,-0.0001564,
	0.0002913,0.0008154,0.0012589,0.0014630,0.0013215,
	0.0008202,0.0000507,-0.0008079,-0.0015377,-0.0019434,
	-0.0019103,-0.0014412,-0.0006566,0.0002388,0.0010162,
	0.0014892,0.0015680,0.0012813,0.0007614,0.0001970,
	-0.0002316,-0.0004128,-0.0003386,-0.0001036,0.0001329,
	0.0002089,0.0000241,-0.0004148,-0.0009865,-0.0014884,
	-0.0017028,-0.0014730,-0.0007676,0.0002944,0.0014645,
	0.0024307,0.0029088,0.0027301,0.0018995,0.0006023,
	-0.0008434,-0.0020756,-0.0027923,-0.0028384,-0.0022480,
	-0.0012285,-0.0000922,0.0008433,0.0013541,0.0013759,
	0.0010148,0.0005013,0.0001008,0.0000139,0.0002974,
	0.0008381,0.0013873,0.0016481,0.0013857,0.0005249,
	-0.0008009,-0.0022692,-0.0034555,-0.0039604,-0.0035413,
	-0.0022027,-0.0002144,0.0019508,0.0037465,0.0047089,
	0.0045927,0.0034451,0.0015917,-0.0004616,-0.0021807,
	-0.0031622,-0.0032487,-0.0025635,-0.0014551,-0.0003677,
	0.0003151,0.0004054,-0.0000257,-0.0006766,-0.0011289,
	-0.0010090,-0.0001413,0.0013604,0.0030996,0.0045058,
	0.0050151,0.0042652,0.0022414,-0.0006794,-0.0038081,
	-0.0063273,-0.0075280,-0.0070257,-0.0048883,-0.0016313,
	0.0019227,0.0048787,0.0065293,0.0065563,0.0051142,
	0.0027637,0.0002808,-0.0016005,-0.0024254,-0.0021589,
	-0.0011796,-0.0001287,0.0003368,-0.0001963,-0.0017160,
	-0.0037507,-0.0055011,-0.0061048,-0.0049487,-0.0019277,
	0.0024461,0.0071303,0.0108241,0.0123497,0.0110281,
	0.0069290,0.0009066,-0.0056017,-0.0109646,-0.0138213,
	-0.0134800,-0.0101259,-0.0047664,0.0010757,0.0058289,
	0.0083564,0.0083005,0.0061722,0.0031516,0.0006560,
	-0.0001855,0.0010254,0.0037532,0.0066543,0.0080120,
	0.0063591,0.0010896,-0.0071501,-0.0164488,-0.0240462,
	-0.0270926,-0.0235498,-0.0129653,0.0031341,0.0213937,
	0.0373397,0.0464959,0.0456227,0.0337361,0.0125992,
	-0.0134918,-0.0386213,-0.0566646,-0.0628429,-0.0550152,
	-0.0343514,-0.0051803,0.0259740,0.0519421,0.0666580,
	0.0666580,0.0519421,0.0259740,-0.0051803,-0.0343514,
	-0.0550152,-0.0628429,-0.0566646,-0.0386213,-0.0134918,
	0.0125992,0.0337361,0.0456227,0.0464959,0.0373397,
	0.0213937,0.0031341,-0.0129653,-0.0235498,-0.0270926,
	-0.0240462,-0.0164488,-0.0071501,0.0010896,0.0063591,
	0.0080120,0.0066543,0.0037532,0.0010254,-0.0001855,
	0.0006560,0.0031516,0.0061722,0.0083005,0.0083564,
	0.0058289,0.0010757,-0.0047664,-0.0101259,-0.0134800,
	-0.0138213,-0.0109646,-0.0056017,0.0009066,0.0069290,
	0.0110281,0.0123497,0.0108241,0.0071303,0.0024461,
	-0.0019277,-0.0049487,-0.0061048,-0.0055011,-0.0037507,
	-0.0017160,-0.0001963,0.0003368,-0.0001287,-0.0011796,
	-0.0021589,-0.0024254,-0.0016005,0.0002808,0.0027637,
	0.0051142,0.0065563,0.0065293,0.0048787,0.0019227,
	-0.0016313,-0.0048883,-0.0070257,-0.0075280,-0.0063273,
	-0.0038081,-0.0006794,0.0022414,0.0042652,0.0050151,
	0.0045058,0.0030996,0.0013604,-0.0001413,-0.0010090,
	-0.0011289,-0.0006766,-0.0000257,0.0004054,0.0003151,
	-0.0003677,-0.0014551,-0.0025635,-0.0032487,-0.0031622,
	-0.0021807,-0.0004616,0.0015917,0.0034451,0.0045927,
	0.0047089,0.0037465,0.0019508,-0.0002144,-0.0022027,
	-0.0035413,-0.0039604,-0.0034555,-0.0022692,-0.0008009,
	0.0005249,0.0013857,0.0016481,0.0013873,0.0008381,
	0.0002974,0.0000139,0.0001008,0.0005013,0.0010148,
	0.0013759,0.0013541,0.0008433,-0.0000922,-0.0012285,
	-0.0022480,-0.0028384,-0.0027923,-0.0020756,-0.0008434,
	0.0006023,0.0018995,0.0027301,0.0029088,0.0024307,
	0.0014645,0.0002944,-0.0007676,-0.0014730,-0.0017028,
	-0.0014884,-0.0009865,-0.0004148,0.0000241,0.0002089,
	0.0001329,-0.0001036,-0.0003386,-0.0004128,-0.0002316,
	0.0001970,0.0007614,0.0012813,0.0015680,0.0014892,
	0.0010162,0.0002388,-0.0006566,-0.0014412,-0.0019103,
	-0.0019434,-0.0015377,-0.0008079,0.0000507,0.0008202,
	0.0013215,0.0014630,0.0012589,0.0008154,0.0002913,
	-0.0001564,-0.0004187,-0.0004646,-0.0003412,-0.0001489,
	-0.0000005,0.0000214,-0.0001057,-0.0003370,-0.0005774,
	-0.0007163,-0.0006690,-0.0004081,0.0000229,0.0005196,
	0.0009471,0.0011822,0.0011506,0.0008504,0.0003541,
	-0.0002128,-0.0007085,-0.0010154,-0.0010713,-0.0008855,
	-0.0005302,-0.0001153,0.0002447,0.0004675,0.0005211,
	0.0004305,0.0002612,0.0000946,-0.0000028,-0.0000004,
	0.0000868,0.0002072,0.0002915,0.0002802,0.0001455,
	-0.0000946,-0.0003807,-0.0006285,-0.0007570,-0.0007129,
	-0.0004890,-0.0001279,0.0002852,0.0006481,0.0008704,
	0.0008985,0.0007234,0.0003947,-0.0000064,-0.0003830,
	-0.0006480,-0.0007505,-0.0006802,-0.0004722,-0.0001896,
	0.0000927,0.0003104,0.0004232,0.0004238,0.0003329,
	0.0001908,0.0000423,-0.0000752,-0.0001415,-0.0001546,
	-0.0001285,-0.0000843,-0.0000436,-0.0000221,-0.0000293,
	-0.0000714,-0.0001608,-0.0003279,-0.0006407,0.0004589	
};

/*
 * Band pass filter coefficients with (900 Hz - 1200 Hz) bandwidth
 */
const float pm FILTER_4[NUM_TAPS] = {
	-0.0002798,0.0007866,0.0002585,0.0000827,0.0000235,
	0.0000150,0.0000333,0.0000540,0.0000474,-0.0000058,
	-0.0000966,-0.0001821,-0.0002032,-0.0001173,0.0000671,
	0.0002809,0.0004166,0.0003807,0.0001487,-0.0002056,
	-0.0005313,-0.0006640,-0.0005096,-0.0001015,0.0003991,
	0.0007699,0.0008286,0.0005269,-0.0000194,-0.0005741,
	-0.0008915,-0.0008282,-0.0004169,0.0001544,0.0006290,
	0.0008060,0.0006360,0.0002264,-0.0002098,-0.0004774,
	-0.0004898,-0.0002996,-0.0000531,0.0001013,0.0001026,
	0.0000054,-0.0000647,-0.0000036,0.0001945,0.0004184,
	0.0005006,0.0003208,-0.0001036,-0.0006012,-0.0009192,
	-0.0008564,-0.0003814,0.0003282,0.0009576,0.0012048,
	0.0009348,0.0002548,-0.0005314,-0.0010705,-0.0011319,
	-0.0007187,-0.0000545,0.0005374,0.0008048,0.0006851,
	0.0003195,-0.0000481,-0.0002216,-0.0001630,-0.0000011,
	0.0000654,-0.0000868,-0.0004080,-0.0006876,-0.0006757,
	-0.0002498,0.0004762,0.0011741,0.0014571,0.0010928,
	0.0001549,-0.0009758,-0.0017784,-0.0018478,-0.0011087,
	0.0001299,0.0013148,0.0019158,0.0016842,0.0007685,
	-0.0003728,-0.0012122,-0.0014191,-0.0010076,-0.0002935,
	0.0003077,0.0005228,0.0003608,0.0000763,-0.0000200,
	0.0002176,0.0006537,0.0009405,0.0007380,-0.0000424,
	-0.0011235,-0.0019622,-0.0020290,-0.0011195,0.0004792,
	0.0020783,0.0029058,0.0025001,0.0009688,-0.0010333,
	-0.0026075,-0.0030481,-0.0021874,-0.0004697,0.0012847,
	0.0022994,0.0022199,0.0012555,0.0000162,-0.0008520,
	-0.0010218,-0.0006311,-0.0001435,-0.0000170,-0.0003856,
	-0.0009543,-0.0011791,-0.0006334,0.0006609,0.0021516,
	0.0030020,0.0025522,0.0007411,-0.0017597,-0.0038118,
	-0.0043614,-0.0029988,-0.0002208,0.0027625,0.0045991,
	0.0044571,0.0024365,-0.0004883,-0.0029665,-0.0039463,
	-0.0031732,-0.0012532,0.0007354,0.0018531,0.0017877,
	0.0009509,0.0001445,0.0000097,0.0005967,0.0013244,
	0.0013641,0.0002252,-0.0018214,-0.0037665,-0.0043673,
	-0.0028776,0.0004097,0.0041526,0.0065455,0.0062319,
	0.0030671,-0.0016979,-0.0059510,-0.0077296,-0.0061958,
	-0.0020666,0.0027382,0.0060539,0.0065364,0.0042784,
	0.0006624,-0.0024492,-0.0037249,-0.0030157,-0.0012965,
	0.0000416,0.0001123,-0.0008944,-0.0018461,-0.0014441,
	0.0008211,0.0041360,0.0066382,0.0064200,0.0027295,
	-0.0033333,-0.0091061,-0.0115948,-0.0090517,-0.0020607,
	0.0065155,0.0127937,0.0137515,0.0087746,0.0000436,
	-0.0084868,-0.0130225,-0.0117855,-0.0058210,0.0016648,
	0.0070679,0.0082826,0.0056881,0.0016915,-0.0008880,
	-0.0006334,0.0015564,0.0030528,0.0012997,-0.0042201,
	-0.0111055,-0.0149629,-0.0118334,-0.0008690,0.0143259,
	0.0267179,0.0290489,0.0177773,-0.0043074,-0.0284132,
	-0.0432652,-0.0404028,-0.0187503,0.0139456,0.0438629,
	0.0570855,0.0461716,0.0141655,-0.0260627,-0.0571383,
	-0.0649013,-0.0450061,-0.0052630,0.0375169,0.0648613,
	0.0648613,0.0375169,-0.0052630,-0.0450061,-0.0649013,
	-0.0571383,-0.0260627,0.0141655,0.0461716,0.0570855,
	0.0438629,0.0139456,-0.0187503,-0.0404028,-0.0432652,
	-0.0284132,-0.0043074,0.0177773,0.0290489,0.0267179,
	0.0143259,-0.0008690,-0.0118334,-0.0149629,-0.0111055,
	-0.0042201,0.0012997,0.0030528,0.0015564,-0.0006334,
	-0.0008880,0.0016915,0.0056881,0.0082826,0.0070679,
	0.0016648,-0.0058210,-0.0117855,-0.0130225,-0.0084868,
	0.0000436,0.0087746,0.0137515,0.0127937,0.0065155,
	-0.0020607,-0.0090517,-0.0115948,-0.0091061,-0.0033333,
	0.0027295,0.0064200,0.0066382,0.0041360,0.0008211,
	-0.0014441,-0.0018461,-0.0008944,0.0001123,0.0000416,
	-0.0012965,-0.0030157,-0.0037249,-0.0024492,0.0006624,
	0.0042784,0.0065364,0.0060539,0.0027382,-0.0020666,
	-0.0061958,-0.0077296,-0.0059510,-0.0016979,0.0030671,
	0.0062319,0.0065455,0.0041526,0.0004097,-0.0028776,
	-0.0043673,-0.0037665,-0.0018214,0.0002252,0.0013641,
	0.0013244,0.0005967,0.0000097,0.0001445,0.0009509,
	0.0017877,0.0018531,0.0007354,-0.0012532,-0.0031732,
	-0.0039463,-0.0029665,-0.0004883,0.0024365,0.0044571,
	0.0045991,0.0027625,-0.0002208,-0.0029988,-0.0043614,
	-0.0038118,-0.0017597,0.0007411,0.0025522,0.0030020,
	0.0021516,0.0006609,-0.0006334,-0.0011791,-0.0009543,
	-0.0003856,-0.0000170,-0.0001435,-0.0006311,-0.0010218,
	-0.0008520,0.0000162,0.0012555,0.0022199,0.0022994,
	0.0012847,-0.0004697,-0.0021874,-0.0030481,-0.0026075,
	-0.0010333,0.0009688,0.0025001,0.0029058,0.0020783,
	0.0004792,-0.0011195,-0.0020290,-0.0019622,-0.0011235,
	-0.0000424,0.0007380,0.0009405,0.0006537,0.0002176,
	-0.0000200,0.0000763,0.0003608,0.0005228,0.0003077,
	-0.0002935,-0.0010076,-0.0014191,-0.0012122,-0.0003728,
	0.0007685,0.0016842,0.0019158,0.0013148,0.0001299,
	-0.0011087,-0.0018478,-0.0017784,-0.0009758,0.0001549,
	0.0010928,0.0014571,0.0011741,0.0004762,-0.0002498,
	-0.0006757,-0.0006876,-0.0004080,-0.0000868,0.0000654,
	-0.0000011,-0.0001630,-0.0002216,-0.0000481,0.0003195,
	0.0006851,0.0008048,0.0005374,-0.0000545,-0.0007187,
	-0.0011319,-0.0010705,-0.0005314,0.0002548,0.0009348,
	0.0012048,0.0009576,0.0003282,-0.0003814,-0.0008564,
	-0.0009192,-0.0006012,-0.0001036,0.0003208,0.0005006,
	0.0004184,0.0001945,-0.0000036,-0.0000647,0.0000054,
	0.0001026,0.0001013,-0.0000531,-0.0002996,-0.0004898,
	-0.0004774,-0.0002098,0.0002264,0.0006360,0.0008060,
	0.0006290,0.0001544,-0.0004169,-0.0008282,-0.0008915,
	-0.0005741,-0.0000194,0.0005269,0.0008286,0.0007699,
	0.0003991,-0.0001015,-0.0005096,-0.0006640,-0.0005313,
	-0.0002056,0.0001487,0.0003807,0.0004166,0.0002809,
	0.0000671,-0.0001173,-0.0002032,-0.0001821,-0.0000966,
	-0.0000058,0.0000474,0.0000540,0.0000333,0.0000150,
	0.0000235,0.0000827,0.0002585,0.0007866,-0.0002798
};

/*
 * Band pass filter coefficients with (1200 Hz - 1500 Hz) bandwidth
 */
const float pm FILTER_5[NUM_TAPS] = {
	0.0001154,-0.0008631,-0.0001296,-0.0000156,0.0000051,
	-0.0000054,-0.0000160,0.0000063,0.0000623,0.0001020,
	0.0000615,-0.0000697,-0.0002101,-0.0002306,-0.0000614,
	0.0002198,0.0004123,0.0003315,-0.0000280,-0.0004457,
	-0.0006083,-0.0003448,0.0002141,0.0006913,0.0007234,
	0.0002397,-0.0004537,-0.0008724,-0.0007015,-0.0000388,
	0.0006612,0.0009093,0.0005377,-0.0001815,-0.0007394,
	-0.0007670,-0.0002905,0.0003135,0.0006231,0.0004786,
	0.0000690,-0.0002656,-0.0003156,-0.0001457,0.0000116,
	0.0000067,-0.0000975,-0.0000986,0.0001136,0.0004042,
	0.0004775,0.0001513,-0.0004270,-0.0008325,-0.0006847,
	0.0000138,0.0008123,0.0011171,0.0006474,-0.0003140,
	-0.0010982,-0.0011385,-0.0004005,0.0005876,0.0011305,
	0.0008835,0.0000866,-0.0006533,-0.0008444,-0.0004614,
	0.0001003,0.0003985,0.0003084,0.0000693,-0.0000012,
	0.0001586,0.0002907,0.0000931,-0.0004244,-0.0008570,
	-0.0007211,0.0000807,0.0010600,0.0014496,0.0008106,
	-0.0005370,-0.0016655,-0.0017067,-0.0005373,0.0010668,
	0.0019716,0.0015270,0.0000576,-0.0013825,-0.0018031,
	-0.0010008,0.0003446,0.0012465,0.0011739,0.0003879,
	-0.0003807,-0.0005904,-0.0003025,-0.0000139,-0.0001030,
	-0.0004277,-0.0004659,0.0001198,0.0010312,0.0014683,
	0.0008093,-0.0007334,-0.0021024,-0.0021469,-0.0005856,
	0.0016249,0.0028990,0.0022076,-0.0000778,-0.0023824,
	-0.0030686,-0.0016587,0.0008108,0.0025821,0.0024948,
	0.0007975,-0.0011482,-0.0019878,-0.0013793,-0.0001019,
	0.0007420,0.0006857,0.0001827,0.0000201,0.0004555,
	0.0009145,0.0005649,-0.0007501,-0.0021268,-0.0022328,
	-0.0005086,0.0021090,0.0036877,0.0027688,-0.0003392,
	-0.0035664,-0.0045260,-0.0023462,0.0015574,0.0044470,
	0.0042805,0.0012326,-0.0024639,-0.0042275,-0.0030377,
	-0.0000644,0.0024058,0.0027996,0.0013421,-0.0004124,
	-0.0010737,-0.0005725,0.0000079,-0.0003131,-0.0013029,
	-0.0016536,-0.0003090,0.0022374,0.0039883,0.0030131,
	-0.0006969,-0.0047665,-0.0060113,-0.0029723,0.0026363,
	0.0069066,0.0065829,0.0016239,-0.0045866,-0.0076652,
	-0.0054825,0.0002718,0.0054276,0.0065082,0.0032343,
	-0.0015506,-0.0043343,-0.0036675,-0.0009620,0.0011345,
	0.0012041,0.0001407,-0.0000550,0.0014457,0.0031828,
	0.0026610,-0.0010377,-0.0057261,-0.0073820,-0.0035073,
	0.0041904,0.0103584,0.0098377,0.0019546,-0.0083028,
	-0.0135623,-0.0095735,0.0013011,0.0115190,0.0138438,
	0.0067417,-0.0045126,-0.0119020,-0.0107071,-0.0027759,
	0.0054367,0.0082547,0.0050619,0.0000289,-0.0021872,
	-0.0007661,0.0008837,-0.0009775,-0.0059237,-0.0088014,
	-0.0042660,0.0072596,0.0178929,0.0175116,0.0025519,
	-0.0189033,-0.0311143,-0.0221522,0.0054760,0.0340637,
	0.0420732,0.0203358,-0.0191073,-0.0494147,-0.0474398,
	-0.0114124,0.0357244,0.0611014,0.0451904,-0.0031539,
	-0.0514888,-0.0659455,-0.0353556,0.0201100,0.0624866,
	0.0624866,0.0201100,-0.0353556,-0.0659455,-0.0514888,
	-0.0031539,0.0451904,0.0611014,0.0357244,-0.0114124,
	-0.0474398,-0.0494147,-0.0191073,0.0203358,0.0420732,
	0.0340637,0.0054760,-0.0221522,-0.0311143,-0.0189033,
	0.0025519,0.0175116,0.0178929,0.0072596,-0.0042660,
	-0.0088014,-0.0059237,-0.0009775,0.0008837,-0.0007661,
	-0.0021872,0.0000289,0.0050619,0.0082547,0.0054367,
	-0.0027759,-0.0107071,-0.0119020,-0.0045126,0.0067417,
	0.0138438,0.0115190,0.0013011,-0.0095735,-0.0135623,
	-0.0083028,0.0019546,0.0098377,0.0103584,0.0041904,
	-0.0035073,-0.0073820,-0.0057261,-0.0010377,0.0026610,
	0.0031828,0.0014457,-0.0000550,0.0001407,0.0012041,
	0.0011345,-0.0009620,-0.0036675,-0.0043343,-0.0015506,
	0.0032343,0.0065082,0.0054276,0.0002718,-0.0054825,
	-0.0076652,-0.0045866,0.0016239,0.0065829,0.0069066,
	0.0026363,-0.0029723,-0.0060113,-0.0047665,-0.0006969,
	0.0030131,0.0039883,0.0022374,-0.0003090,-0.0016536,
	-0.0013029,-0.0003131,0.0000079,-0.0005725,-0.0010737,
	-0.0004124,0.0013421,0.0027996,0.0024058,-0.0000644,
	-0.0030377,-0.0042275,-0.0024639,0.0012326,0.0042805,
	0.0044470,0.0015574,-0.0023462,-0.0045260,-0.0035664,
	-0.0003392,0.0027688,0.0036877,0.0021090,-0.0005086,
	-0.0022328,-0.0021268,-0.0007501,0.0005649,0.0009145,
	0.0004555,0.0000201,0.0001827,0.0006857,0.0007420,
	-0.0001019,-0.0013793,-0.0019878,-0.0011482,0.0007975,
	0.0024948,0.0025821,0.0008108,-0.0016587,-0.0030686,
	-0.0023824,-0.0000778,0.0022076,0.0028990,0.0016249,
	-0.0005856,-0.0021469,-0.0021024,-0.0007334,0.0008093,
	0.0014683,0.0010312,0.0001198,-0.0004659,-0.0004277,
	-0.0001030,-0.0000139,-0.0003025,-0.0005904,-0.0003807,
	0.0003879,0.0011739,0.0012465,0.0003446,-0.0010008,
	-0.0018031,-0.0013825,0.0000576,0.0015270,0.0019716,
	0.0010668,-0.0005373,-0.0017067,-0.0016655,-0.0005370,
	0.0008106,0.0014496,0.0010600,0.0000807,-0.0007211,
	-0.0008570,-0.0004244,0.0000931,0.0002907,0.0001586,
	-0.0000012,0.0000693,0.0003084,0.0003985,0.0001003,
	-0.0004614,-0.0008444,-0.0006533,0.0000866,0.0008835,
	0.0011305,0.0005876,-0.0004005,-0.0011385,-0.0010982,
	-0.0003140,0.0006474,0.0011171,0.0008123,0.0000138,
	-0.0006847,-0.0008325,-0.0004270,0.0001513,0.0004775,
	0.0004042,0.0001136,-0.0000986,-0.0000975,0.0000067,
	0.0000116,-0.0001457,-0.0003156,-0.0002656,0.0000690,
	0.0004786,0.0006231,0.0003135,-0.0002905,-0.0007670,
	-0.0007394,-0.0001815,0.0005377,0.0009093,0.0006612,
	-0.0000388,-0.0007015,-0.0008724,-0.0004537,0.0002397,
	0.0007234,0.0006913,0.0002141,-0.0003448,-0.0006083,
	-0.0004457,-0.0000280,0.0003315,0.0004123,0.0002198,
	-0.0000614,-0.0002306,-0.0002101,-0.0000697,0.0000615,
	0.0001020,0.0000623,0.0000063,-0.0000160,-0.0000054,
	0.0000051,-0.0000156,-0.0001296,-0.0008631,0.0001154
};