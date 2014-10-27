#ifndef FIBERMAP_H
#define FIBERMAP_H

// look-up table of detector fibers to pade channels
// naming conventions
//   PADE channels - boardID*100+channel_number
//   FIBER numbers - moduleID*100+fiber_number (+/- for front/back)
const int NPADECHANNELS=128;
// board/chan	module/fiber(Down Str)	board/chan   module/fiber(UP Stream)
// ChannelID    FiberID                 ChannelID    FiberID                            
const int FIBERMAP_OCTOBER14[]={
//	Cookie A			Cookie B		
	11206,	101,		11525,		-101,
	11204,	102,		11527,		-102,
	11203,	103,		11528,		-103,
	11201,	104,		11530,		-104,
	11214,	201,		11517,		-201,
	11212,	202,		11519,		-202,
	11211,	203,		11520,		-203,
	11209,	204,		11522,		-204,
	11217,	301,		11514,		-301,
	11219,	302,		11512,		-302,
	11220,	303,		11511,		-303,
	11222,	304,		11509,		-304,
	11225,	401,		11506,		-401,
	11227,	402,		11504,		-402,
	11228,	403,		11501,		-403,
	11230,	404,		11503,		-404,
	11200,	501,		11531,		-501,
	11202,	502,		11529,		-502,
	11205,	503,		11526,		-503,
	11207,	504,		11524,		-504,
	11208,	601,		11523,		-601,
	11210,	602,		11521,		-602,
	11213,	603,		11518,		-603,
	11215,	604,		11516,		-604,
	11223,	701,		11508,		-701,
	11221,	702,		11510,		-702,
	11218,	703,		11513,		-703,
	11216,	704,		11515,		-704,
	11231,	801,		11500,		-801,
	11229,	802,		11502,		-802,
	11226,	803,		11507,		-803,
	11224,	804,		11505,		-804,
	1606,	901,		11625,		-901,
	1604,	902,		11627,		-902,
	1603,	903,		11628,		-903,
	1601,	904,		11630,		-904,
	1614,	1001,		11617,		-1001,
	1612,	1002,		11619,		-1002,
	1611,	1003,		11620,		-1003,
	1609,	1004,		11622,		-1004,
	1617,	1101,		11614,		-1101,
	1619,	1102,		11612,		-1102,
	1620,	1103,		11611,		-1103,
	1622,	1104,		11609,		-1104,
	1625,	1201,		11606,		-1201,
	1627,	1202,		11604,		-1202,
	1628,	1203,		11603,		-1203,
	1630,	1204,		11601,		-1204,
	1600,	1301,		11631,		-1301,
	1602,	1302,		11629,		-1302,
	1605,	1303,		11624,		-1303,
	1607,	1304,		11626,		-1304,
	1608,	1401,		11623,		-1401,
	1610,	1402,		11621,		-1402,
	1613,	1403,		11618,		-1403,
	1615,	1404,		11616,		-1404,
	1623,	1501,		11608,		-1501,
	1621,	1502,		11610,		-1502,
	1618,	1503,		11615,		-1503,
	1616,	1504,		11613,		-1504,
	1631,	1601,		11600,		-1601,
	1629,	1602,		11602,		-1602,
	1626,	1603,		11605,		-1603,
	1624,	1604,		11607,		-1604
};


const int FIBERMAP_APRIL14[]={
//	Cookie A			Cookie B		
	11506,	101,		11225,		-101,
	11504,	102,		11227,		-102,
	11503,	103,		11228,		-103,
	11501,	104,		11230,		-104,
	11514,	201,		11217,		-201,
	11512,	202,		11219,		-202,
	11511,	203,		11220,		-203,
	11509,	204,		11222,		-204,
	11517,	301,		11214,		-301,
	11519,	302,		11212,		-302,
	11520,	303,		11211,		-303,
	11522,	304,		11209,		-304,
	11525,	401,		11206,		-401,
	11527,	402,		11204,		-402,
	11528,	403,		11201,		-403,
	11530,	404,		11203,		-404,
	11500,	501,		11231,		-501,
	11502,	502,		11229,		-502,
	11505,	503,		11226,		-503,
	11507,	504,		11224,		-504,
	11508,	601,		11223,		-601,
	11510,	602,		11221,		-602,
	11513,	603,		11218,		-603,
	11515,	604,		11216,		-604,
	11523,	701,		11208,		-701,
	11521,	702,		11210,		-702,
	11518,	703,		11213,		-703,
	11516,	704,		11215,		-704,
	11531,	801,		11200,		-801,
	11529,	802,		11202,		-802,
	11526,	803,		11207,		-803,
	11524,	804,		11205,		-804,
	11306,	901,		11625,		-901,
	11304,	902,		11627,		-902,
	11303,	903,		11628,		-903,
	11301,	904,		11630,		-904,
	11314,	1001,		11617,		-1001,
	11312,	1002,		11619,		-1002,
	11311,	1003,		11620,		-1003,
	11309,	1004,		11622,		-1004,
	11317,	1101,		11614,		-1101,
	11319,	1102,		11612,		-1102,
	11320,	1103,		11611,		-1103,
	11322,	1104,		11609,		-1104,
	11325,	1201,		11606,		-1201,
	11327,	1202,		11604,		-1202,
	11328,	1203,		11603,		-1203,
	11330,	1204,		11601,		-1204,
	11300,	1301,		11631,		-1301,
	11302,	1302,		11629,		-1302,
	11305,	1303,		11624,		-1303,
	11307,	1304,		11626,		-1304,
	11308,	1401,		11623,		-1401,
	11310,	1402,		11621,		-1402,
	11313,	1403,		11618,		-1403,
	11315,	1404,		11616,		-1404,
	11323,	1501,		11608,		-1501,
	11321,	1502,		11610,		-1502,
	11318,	1503,		11615,		-1503,
	11316,	1504,		11613,		-1504,
	11331,	1601,		11600,		-1601,
	11329,	1602,		11602,		-1602,
	11326,	1603,		11605,		-1603,
	11324,	1604,		11607,		-1604
};


const int FIBERMAP_JULY14[]={
//	Cookie A			Cookie B		
	11206,	101,		11525,		-101,
	11204,	102,		11527,		-102,
	11203,	103,		11528,		-103,
	11201,	104,		11530,		-104,
	11214,	201,		11517,		-201,
	11212,	202,		11519,		-202,
	11211,	203,		11520,		-203,
	11209,	204,		11522,		-204,
	11217,	301,		11514,		-301,
	11219,	302,		11512,		-302,
	11220,	303,		11511,		-303,
	11222,	304,		11509,		-304,
	11225,	401,		11506,		-401,
	11227,	402,		11504,		-402,
	11228,	403,		11501,		-403,
	11230,	404,		11503,		-404,
	11200,	501,		11531,		-501,
	11202,	502,		11529,		-502,
	11205,	503,		11526,		-503,
	11207,	504,		11524,		-504,
	11208,	601,		11523,		-601,
	11210,	602,		11521,		-602,
	11213,	603,		11518,		-603,
	11215,	604,		11516,		-604,
	11223,	701,		11508,		-701,
	11221,	702,		11510,		-702,
	11218,	703,		11513,		-703,
	11216,	704,		11515,		-704,
	11231,	801,		11500,		-801,
	11229,	802,		11502,		-802,
	11226,	803,		11507,		-803,
	11224,	804,		11505,		-804,
	11706,	901,		11625,		-901,
	11704,	902,		11627,		-902,
	11703,	903,		11628,		-903,
	11701,	904,		11630,		-904,
	11714,	1001,		11617,		-1001,
	11712,	1002,		11619,		-1002,
	11711,	1003,		11620,		-1003,
	11709,	1004,		11622,		-1004,
	11717,	1101,		11614,		-1101,
	11719,	1102,		11612,		-1102,
	11720,	1103,		11611,		-1103,
	11722,	1104,		11609,		-1104,
	11725,	1201,		11606,		-1201,
	11727,	1202,		11604,		-1202,
	11728,	1203,		11603,		-1203,
	11730,	1204,		11601,		-1204,
	11700,	1301,		11631,		-1301,
	11702,	1302,		11629,		-1302,
	11705,	1303,		11624,		-1303,
	11707,	1304,		11626,		-1304,
	11708,	1401,		11623,		-1401,
	11710,	1402,		11621,		-1402,
	11713,	1403,		11618,		-1403,
	11715,	1404,		11616,		-1404,
	11723,	1501,		11608,		-1501,
	11721,	1502,		11610,		-1502,
	11718,	1503,		11615,		-1503,
	11716,	1504,		11613,		-1504,
	11731,	1601,		11600,		-1601,
	11729,	1602,		11602,		-1602,
	11726,	1603,		11605,		-1603,
	11724,	1604,		11607,		-1604
};

#endif