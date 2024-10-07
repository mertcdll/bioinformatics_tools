import pandas as pd
import json

def getBaseINT(bases:str) -> str:

    ints = {
        "A":"1",
        "T":"2",
        "C":"3",
        "G":"4"
    }
    res = ""
    for base in bases:
        try:
            res += ints[base.upper()]
        except Exception as e:
            pass # handle the exception later
    return res

def getChromINT(chrom:str) -> str:

    ints = {"chr1": "1111", "chr2": "1112", "chr3": "1113", "chr4": "1114", "chr5": "1115", "chr6": "1116", "chr7": "1117", "chr8": "1118", "chr9": "1119", "chr10": "1120", "chr11": "1121", "chr12": "1122", "chr13": "1123", "chr14": "1124", "chr15": "1125", "chr16": "1126", "chr17": "1127", "chr18": "1128", "chr19": "1129", "chr20": "1130", "chr21": "1131", "chr22": "1132", "chrX": "1133", "chrY": "1134", "chrM": "1135", "chr1_GL383518v1_alt": "1136", "chr1_GL383519v1_alt": "1137", "chr1_GL383520v2_alt": "1138", "chr1_KI270706v1_random": "1139", "chr1_KI270707v1_random": "1140", "chr1_KI270708v1_random": "1141", "chr1_KI270709v1_random": "1142", "chr1_KI270710v1_random": "1143", "chr1_KI270711v1_random": "1144", "chr1_KI270712v1_random": "1145", "chr1_KI270713v1_random": "1146", "chr1_KI270714v1_random": "1147", "chr1_KI270759v1_alt": "1148", "chr1_KI270760v1_alt": "1149", "chr1_KI270761v1_alt": "1150", "chr1_KI270762v1_alt": "1151", "chr1_KI270763v1_alt": "1152", "chr1_KI270764v1_alt": "1153", "chr1_KI270765v1_alt": "1154", "chr1_KI270766v1_alt": "1155", "chr1_KI270892v1_alt": "1156", "chr1_KN196472v1_fix": "1157", "chr1_KN196473v1_fix": "1158", "chr1_KN196474v1_fix": "1159", "chr1_KN538360v1_fix": "1160", "chr1_KN538361v1_fix": "1161", "chr1_KQ031383v1_fix": "1162", "chr1_KQ458382v1_alt": "1163", "chr1_KQ458383v1_alt": "1164", "chr1_KQ458384v1_alt": "1165", "chr1_KQ983255v1_alt": "1166", "chr1_KV880763v1_alt": "1167", "chr1_KZ208904v1_alt": "1168", "chr1_KZ208905v1_alt": "1169", "chr1_KZ208906v1_fix": "1170", "chr1_KZ559100v1_fix": "1171", "chr1_MU273330v1_alt": "1172", "chr1_MU273331v1_alt": "1173", "chr1_MU273332v1_alt": "1174", "chr1_MU273333v1_fix": "1175", "chr1_MU273334v1_fix": "1176", "chr1_MU273335v1_fix": "1177", "chr1_MU273336v1_fix": "1178", "chr2_GL383521v1_alt": "1179", "chr2_GL383522v1_alt": "1180", "chr2_GL582966v2_alt": "1181", "chr2_KI270715v1_random": "1182", "chr2_KI270716v1_random": "1183", "chr2_KI270767v1_alt": "1184", "chr2_KI270768v1_alt": "1185", "chr2_KI270769v1_alt": "1186", "chr2_KI270770v1_alt": "1187", "chr2_KI270771v1_alt": "1188", "chr2_KI270772v1_alt": "1189", "chr2_KI270773v1_alt": "1190", "chr2_KI270774v1_alt": "1191", "chr2_KI270775v1_alt": "1192", "chr2_KI270776v1_alt": "1193", "chr2_KI270893v1_alt": "1194", "chr2_KI270894v1_alt": "1195", "chr2_KN538362v1_fix": "1196", "chr2_KN538363v1_fix": "1197", "chr2_KQ031384v1_fix": "1198", "chr2_KQ983256v1_alt": "1199", "chr2_KZ208907v1_alt": "1200", "chr2_KZ208908v1_alt": "1201", "chr2_ML143341v1_fix": "1202", "chr2_ML143342v1_fix": "1203", "chr2_MU273337v1_alt": "1204", "chr2_MU273338v1_alt": "1205", "chr2_MU273339v1_alt": "1206", "chr2_MU273340v1_alt": "1207", "chr2_MU273341v1_fix": "1208", "chr2_MU273342v1_fix": "1209", "chr2_MU273343v1_fix": "1210", "chr2_MU273344v1_fix": "1211", "chr2_MU273345v1_fix": "1212", "chr3_GL000221v1_random": "1213", "chr3_GL383526v1_alt": "1214", "chr3_JH636055v2_alt": "1215", "chr3_KI270777v1_alt": "1216", "chr3_KI270778v1_alt": "1217", "chr3_KI270779v1_alt": "1218", "chr3_KI270780v1_alt": "1219", "chr3_KI270781v1_alt": "1220", "chr3_KI270782v1_alt": "1221", "chr3_KI270783v1_alt": "1222", "chr3_KI270784v1_alt": "1223", "chr3_KI270895v1_alt": "1224", "chr3_KI270924v1_alt": "1225", "chr3_KI270934v1_alt": "1226", "chr3_KI270935v1_alt": "1227", "chr3_KI270936v1_alt": "1228", "chr3_KI270937v1_alt": "1229", "chr3_KN196475v1_fix": "1230", "chr3_KN196476v1_fix": "1231", "chr3_KN538364v1_fix": "1232", "chr3_KQ031385v1_fix": "1233", "chr3_KQ031386v1_fix": "1234", "chr3_KV766192v1_fix": "1235", "chr3_KZ208909v1_alt": "1236", "chr3_KZ559101v1_alt": "1237", "chr3_KZ559102v1_alt": "1238", "chr3_KZ559103v1_alt": "1239", "chr3_KZ559104v1_fix": "1240", "chr3_KZ559105v1_alt": "1241", "chr3_ML143343v1_alt": "1242", "chr3_MU273346v1_fix": "1243", "chr3_MU273347v1_fix": "1244", "chr3_MU273348v1_fix": "1245", "chr4_GL000008v2_random": "1246", "chr4_GL000257v2_alt": "1247", "chr4_GL383527v1_alt": "1248", "chr4_GL383528v1_alt": "1249", "chr4_KI270785v1_alt": "1250", "chr4_KI270786v1_alt": "1251", "chr4_KI270787v1_alt": "1252", "chr4_KI270788v1_alt": "1253", "chr4_KI270789v1_alt": "1254", "chr4_KI270790v1_alt": "1255", "chr4_KI270896v1_alt": "1256", "chr4_KI270925v1_alt": "1257", "chr4_KQ090013v1_alt": "1258", "chr4_KQ090014v1_alt": "1259", "chr4_KQ090015v1_alt": "1260", "chr4_KQ983257v1_fix": "1261", "chr4_KQ983258v1_alt": "1262", "chr4_KV766193v1_alt": "1263", "chr4_ML143344v1_fix": "1264", "chr4_ML143345v1_fix": "1265", "chr4_ML143346v1_fix": "1266", "chr4_ML143347v1_fix": "1267", "chr4_ML143348v1_fix": "1268", "chr4_ML143349v1_fix": "1269", "chr4_MU273349v1_alt": "1270", "chr4_MU273350v1_fix": "1271", "chr4_MU273351v1_fix": "1272", "chr5_GL000208v1_random": "1273", "chr5_GL339449v2_alt": "1274", "chr5_GL383530v1_alt": "1275", "chr5_GL383531v1_alt": "1276", "chr5_GL383532v1_alt": "1277", "chr5_GL949742v1_alt": "1278", "chr5_KI270791v1_alt": "1279", "chr5_KI270792v1_alt": "1280", "chr5_KI270793v1_alt": "1281", "chr5_KI270794v1_alt": "1282", "chr5_KI270795v1_alt": "1283", "chr5_KI270796v1_alt": "1284", "chr5_KI270897v1_alt": "1285", "chr5_KI270898v1_alt": "1286", "chr5_KN196477v1_alt": "1287", "chr5_KV575243v1_alt": "1288", "chr5_KV575244v1_fix": "1289", "chr5_KZ208910v1_alt": "1290", "chr5_ML143350v1_fix": "1291", "chr5_MU273352v1_fix": "1292", "chr5_MU273353v1_fix": "1293", "chr5_MU273354v1_fix": "1294", "chr5_MU273355v1_fix": "1295", "chr5_MU273356v1_alt": "1296", "chr6_GL000250v2_alt": "1297", "chr6_GL000251v2_alt": "1298", "chr6_GL000252v2_alt": "1299", "chr6_GL000253v2_alt": "1300", "chr6_GL000254v2_alt": "1301", "chr6_GL000255v2_alt": "1302", "chr6_GL000256v2_alt": "1303", "chr6_GL383533v1_alt": "1304", "chr6_KB021644v2_alt": "1305", "chr6_KI270758v1_alt": "1306", "chr6_KI270797v1_alt": "1307", "chr6_KI270798v1_alt": "1308", "chr6_KI270799v1_alt": "1309", "chr6_KI270800v1_alt": "1310", "chr6_KI270801v1_alt": "1311", "chr6_KI270802v1_alt": "1312", "chr6_KN196478v1_fix": "1313", "chr6_KQ031387v1_fix": "1314", "chr6_KQ090016v1_fix": "1315", "chr6_KQ090017v1_alt": "1316", "chr6_KV766194v1_fix": "1317", "chr6_KZ208911v1_fix": "1318", "chr6_ML143351v1_fix": "1319", "chr6_MU273357v1_alt": "1320", "chr7_GL383534v2_alt": "1321", "chr7_KI270803v1_alt": "1322", "chr7_KI270804v1_alt": "1323", "chr7_KI270805v1_alt": "1324", "chr7_KI270806v1_alt": "1325", "chr7_KI270807v1_alt": "1326", "chr7_KI270808v1_alt": "1327", "chr7_KI270809v1_alt": "1328", "chr7_KI270899v1_alt": "1329", "chr7_KQ031388v1_fix": "1330", "chr7_KV880764v1_fix": "1331", "chr7_KV880765v1_fix": "1332", "chr7_KZ208912v1_fix": "1333", "chr7_KZ208913v1_alt": "1334", "chr7_KZ559106v1_alt": "1335", "chr7_ML143352v1_fix": "1336", "chr7_MU273358v1_alt": "1337", "chr8_KI270810v1_alt": "1338", "chr8_KI270811v1_alt": "1339", "chr8_KI270812v1_alt": "1340", "chr8_KI270813v1_alt": "1341", "chr8_KI270814v1_alt": "1342", "chr8_KI270815v1_alt": "1343", "chr8_KI270816v1_alt": "1344", "chr8_KI270817v1_alt": "1345", "chr8_KI270818v1_alt": "1346", "chr8_KI270819v1_alt": "1347", "chr8_KI270820v1_alt": "1348", "chr8_KI270821v1_alt": "1349", "chr8_KI270822v1_alt": "1350", "chr8_KI270900v1_alt": "1351", "chr8_KI270901v1_alt": "1352", "chr8_KI270926v1_alt": "1353", "chr8_KV880766v1_fix": "1354", "chr8_KV880767v1_fix": "1355", "chr8_KZ208914v1_fix": "1356", "chr8_KZ208915v1_fix": "1357", "chr8_KZ559107v1_alt": "1358", "chr8_MU273359v1_fix": "1359", "chr8_MU273360v1_fix": "1360", "chr8_MU273361v1_fix": "1361", "chr8_MU273362v1_fix": "1362", "chr8_MU273363v1_fix": "1363", "chr9_GL383539v1_alt": "1364", "chr9_GL383540v1_alt": "1365", "chr9_GL383541v1_alt": "1366", "chr9_GL383542v1_alt": "1367", "chr9_KI270717v1_random": "1368", "chr9_KI270718v1_random": "1369", "chr9_KI270719v1_random": "1370", "chr9_KI270720v1_random": "1371", "chr9_KI270823v1_alt": "1372", "chr9_KN196479v1_fix": "1373", "chr9_KQ090018v1_alt": "1374", "chr9_KQ090019v1_alt": "1375", "chr9_ML143353v1_fix": "1376", "chr9_MU273364v1_fix": "1377", "chr9_MU273365v1_fix": "1378", "chr9_MU273366v1_fix": "1379", "chr10_GL383545v1_alt": "1380", "chr10_GL383546v1_alt": "1381", "chr10_KI270824v1_alt": "1382", "chr10_KI270825v1_alt": "1383", "chr10_KN196480v1_fix": "1384", "chr10_KN538365v1_fix": "1385", "chr10_KN538366v1_fix": "1386", "chr10_KN538367v1_fix": "1387", "chr10_KQ090020v1_alt": "1388", "chr10_KQ090021v1_fix": "1389", "chr10_ML143354v1_fix": "1390", "chr10_ML143355v1_fix": "1391", "chr10_MU273367v1_fix": "1392", "chr11_GL383547v1_alt": "1393", "chr11_JH159136v1_alt": "1394", "chr11_JH159137v1_alt": "1395", "chr11_KI270721v1_random": "1396", "chr11_KI270826v1_alt": "1397", "chr11_KI270827v1_alt": "1398", "chr11_KI270829v1_alt": "1399", "chr11_KI270830v1_alt": "1400", "chr11_KI270831v1_alt": "1401", "chr11_KI270832v1_alt": "1402", "chr11_KI270902v1_alt": "1403", "chr11_KI270903v1_alt": "1404", "chr11_KI270927v1_alt": "1405", "chr11_KN196481v1_fix": "1406", "chr11_KN538368v1_alt": "1407", "chr11_KQ090022v1_fix": "1408", "chr11_KQ759759v1_fix": "1409", "chr11_KQ759759v2_fix": "1410", "chr11_KV766195v1_fix": "1411", "chr11_KZ559108v1_fix": "1412", "chr11_KZ559109v1_fix": "1413", "chr11_KZ559110v1_alt": "1414", "chr11_KZ559111v1_alt": "1415", "chr11_ML143356v1_fix": "1416", "chr11_ML143357v1_fix": "1417", "chr11_ML143358v1_fix": "1418", "chr11_ML143359v1_fix": "1419", "chr11_ML143360v1_fix": "1420", "chr11_MU273368v1_alt": "1421", "chr11_MU273369v1_fix": "1422", "chr11_MU273370v1_fix": "1423", "chr11_MU273371v1_fix": "1424", "chr12_GL383549v1_alt": "1425", "chr12_GL383550v2_alt": "1426", "chr12_GL383551v1_alt": "1427", "chr12_GL383552v1_alt": "1428", "chr12_GL383553v2_alt": "1429", "chr12_GL877875v1_alt": "1430", "chr12_GL877876v1_alt": "1431", "chr12_KI270833v1_alt": "1432", "chr12_KI270834v1_alt": "1433", "chr12_KI270835v1_alt": "1434", "chr12_KI270836v1_alt": "1435", "chr12_KI270837v1_alt": "1436", "chr12_KI270904v1_alt": "1437", "chr12_KN196482v1_fix": "1438", "chr12_KN538369v1_fix": "1439", "chr12_KN538370v1_fix": "1440", "chr12_KQ090023v1_alt": "1441", "chr12_KQ759760v1_fix": "1442", "chr12_KZ208916v1_fix": "1443", "chr12_KZ208917v1_fix": "1444", "chr12_KZ208918v1_alt": "1445", "chr12_KZ559112v1_alt": "1446", "chr12_ML143361v1_fix": "1447", "chr12_ML143362v1_fix": "1448", "chr12_MU273372v1_fix": "1449", "chr13_KI270838v1_alt": "1450", "chr13_KI270839v1_alt": "1451", "chr13_KI270840v1_alt": "1452", "chr13_KI270841v1_alt": "1453", "chr13_KI270842v1_alt": "1454", "chr13_KI270843v1_alt": "1455", "chr13_KN196483v1_fix": "1456", "chr13_KN538371v1_fix": "1457", "chr13_KN538372v1_fix": "1458", "chr13_KN538373v1_fix": "1459", "chr13_KQ090024v1_alt": "1460", "chr13_KQ090025v1_alt": "1461", "chr13_ML143363v1_fix": "1462", "chr13_ML143364v1_fix": "1463", "chr13_ML143365v1_fix": "1464", "chr13_ML143366v1_fix": "1465", "chr14_GL000009v2_random": "1466", "chr14_GL000194v1_random": "1467", "chr14_GL000225v1_random": "1468", "chr14_KI270722v1_random": "1469", "chr14_KI270723v1_random": "1470", "chr14_KI270724v1_random": "1471", "chr14_KI270725v1_random": "1472", "chr14_KI270726v1_random": "1473", "chr14_KI270844v1_alt": "1474", "chr14_KI270845v1_alt": "1475", "chr14_KI270846v1_alt": "1476", "chr14_KI270847v1_alt": "1477", "chr14_KZ208919v1_alt": "1478", "chr14_KZ208920v1_fix": "1479", "chr14_ML143367v1_fix": "1480", "chr14_ML143368v1_alt": "1481", "chr14_MU273373v1_fix": "1482", "chr15_GL383554v1_alt": "1483", "chr15_GL383555v2_alt": "1484", "chr15_KI270727v1_random": "1485", "chr15_KI270848v1_alt": "1486", "chr15_KI270849v1_alt": "1487", "chr15_KI270850v1_alt": "1488", "chr15_KI270851v1_alt": "1489", "chr15_KI270852v1_alt": "1490", "chr15_KI270905v1_alt": "1491", "chr15_KI270906v1_alt": "1492", "chr15_KN538374v1_fix": "1493", "chr15_KQ031389v1_alt": "1494", "chr15_ML143369v1_fix": "1495", "chr15_ML143370v1_fix": "1496", "chr15_ML143371v1_fix": "1497", "chr15_ML143372v1_fix": "1498", "chr15_MU273374v1_fix": "1499", "chr15_MU273375v1_alt": "1500", "chr16_GL383556v1_alt": "1501", "chr16_GL383557v1_alt": "1502", "chr16_KI270728v1_random": "1503", "chr16_KI270853v1_alt": "1504", "chr16_KI270854v1_alt": "1505", "chr16_KI270855v1_alt": "1506", "chr16_KI270856v1_alt": "1507", "chr16_KQ031390v1_alt": "1508", "chr16_KQ090026v1_alt": "1509", "chr16_KQ090027v1_alt": "1510", "chr16_KV880768v1_fix": "1511", "chr16_KZ208921v1_alt": "1512", "chr16_KZ559113v1_fix": "1513", "chr16_ML143373v1_fix": "1514", "chr16_MU273376v1_fix": "1515", "chr16_MU273377v1_fix": "1516", "chr17_GL000205v2_random": "1517", "chr17_GL000258v2_alt": "1518", "chr17_GL383563v3_alt": "1519", "chr17_GL383564v2_alt": "1520", "chr17_GL383565v1_alt": "1521", "chr17_GL383566v1_alt": "1522", "chr17_JH159146v1_alt": "1523", "chr17_JH159147v1_alt": "1524", "chr17_JH159148v1_alt": "1525", "chr17_KI270729v1_random": "1526", "chr17_KI270730v1_random": "1527", "chr17_KI270857v1_alt": "1528", "chr17_KI270858v1_alt": "1529", "chr17_KI270859v1_alt": "1530", "chr17_KI270860v1_alt": "1531", "chr17_KI270861v1_alt": "1532", "chr17_KI270862v1_alt": "1533", "chr17_KI270907v1_alt": "1534", "chr17_KI270908v1_alt": "1535", "chr17_KI270909v1_alt": "1536", "chr17_KI270910v1_alt": "1537", "chr17_KV575245v1_fix": "1538", "chr17_KV766196v1_fix": "1539", "chr17_KV766197v1_alt": "1540", "chr17_KV766198v1_alt": "1541", "chr17_KZ559114v1_alt": "1542", "chr17_ML143374v1_fix": "1543", "chr17_ML143375v1_fix": "1544", "chr17_MU273378v1_alt": "1545", "chr17_MU273379v1_fix": "1546", "chr17_MU273380v1_fix": "1547", "chr17_MU273381v1_fix": "1548", "chr17_MU273382v1_fix": "1549", "chr17_MU273383v1_fix": "1550", "chr18_GL383567v1_alt": "1551", "chr18_GL383568v1_alt": "1552", "chr18_GL383569v1_alt": "1553", "chr18_GL383570v1_alt": "1554", "chr18_GL383571v1_alt": "1555", "chr18_GL383572v1_alt": "1556", "chr18_KI270863v1_alt": "1557", "chr18_KI270864v1_alt": "1558", "chr18_KI270911v1_alt": "1559", "chr18_KI270912v1_alt": "1560", "chr18_KQ090028v1_fix": "1561", "chr18_KQ458385v1_alt": "1562", "chr18_KZ208922v1_fix": "1563", "chr18_KZ559115v1_fix": "1564", "chr18_KZ559116v1_alt": "1565", "chr19_GL000209v2_alt": "1566", "chr19_GL383573v1_alt": "1567", "chr19_GL383574v1_alt": "1568", "chr19_GL383575v2_alt": "1569", "chr19_GL383576v1_alt": "1570", "chr19_GL949746v1_alt": "1571", "chr19_GL949747v2_alt": "1572", "chr19_GL949748v2_alt": "1573", "chr19_GL949749v2_alt": "1574", "chr19_GL949750v2_alt": "1575", "chr19_GL949751v2_alt": "1576", "chr19_GL949752v1_alt": "1577", "chr19_GL949753v2_alt": "1578", "chr19_KI270865v1_alt": "1579", "chr19_KI270866v1_alt": "1580", "chr19_KI270867v1_alt": "1581", "chr19_KI270868v1_alt": "1582", "chr19_KI270882v1_alt": "1583", "chr19_KI270883v1_alt": "1584", "chr19_KI270884v1_alt": "1585", "chr19_KI270885v1_alt": "1586", "chr19_KI270886v1_alt": "1587", "chr19_KI270887v1_alt": "1588", "chr19_KI270888v1_alt": "1589", "chr19_KI270889v1_alt": "1590", "chr19_KI270890v1_alt": "1591", "chr19_KI270891v1_alt": "1592", "chr19_KI270914v1_alt": "1593", "chr19_KI270915v1_alt": "1594", "chr19_KI270916v1_alt": "1595", "chr19_KI270917v1_alt": "1596", "chr19_KI270918v1_alt": "1597", "chr19_KI270919v1_alt": "1598", "chr19_KI270920v1_alt": "1599", "chr19_KI270921v1_alt": "1600", "chr19_KI270922v1_alt": "1601", "chr19_KI270923v1_alt": "1602", "chr19_KI270929v1_alt": "1603", "chr19_KI270930v1_alt": "1604", "chr19_KI270931v1_alt": "1605", "chr19_KI270932v1_alt": "1606", "chr19_KI270933v1_alt": "1607", "chr19_KI270938v1_alt": "1608", "chr19_KN196484v1_fix": "1609", "chr19_KQ458386v1_fix": "1610", "chr19_KV575246v1_alt": "1611", "chr19_KV575247v1_alt": "1612", "chr19_KV575248v1_alt": "1613", "chr19_KV575249v1_alt": "1614", "chr19_KV575250v1_alt": "1615", "chr19_KV575251v1_alt": "1616", "chr19_KV575252v1_alt": "1617", "chr19_KV575253v1_alt": "1618", "chr19_KV575254v1_alt": "1619", "chr19_KV575255v1_alt": "1620", "chr19_KV575256v1_alt": "1621", "chr19_KV575257v1_alt": "1622", "chr19_KV575258v1_alt": "1623", "chr19_KV575259v1_alt": "1624", "chr19_KV575260v1_alt": "1625", "chr19_ML143376v1_fix": "1626", "chr19_MU273384v1_fix": "1627", "chr19_MU273385v1_fix": "1628", "chr19_MU273386v1_fix": "1629", "chr19_MU273387v1_alt": "1630", "chr20_GL383577v2_alt": "1631", "chr20_KI270869v1_alt": "1632", "chr20_KI270870v1_alt": "1633", "chr20_KI270871v1_alt": "1634", "chr20_MU273388v1_fix": "1635", "chr20_MU273389v1_fix": "1636", "chr21_GL383578v2_alt": "1637", "chr21_GL383579v2_alt": "1638", "chr21_GL383580v2_alt": "1639", "chr21_GL383581v2_alt": "1640", "chr21_KI270872v1_alt": "1641", "chr21_KI270873v1_alt": "1642", "chr21_KI270874v1_alt": "1643", "chr21_ML143377v1_fix": "1644", "chr21_MU273390v1_fix": "1645", "chr21_MU273391v1_fix": "1646", "chr21_MU273392v1_fix": "1647", "chr22_GL383582v2_alt": "1648", "chr22_GL383583v2_alt": "1649", "chr22_KB663609v1_alt": "1650", "chr22_KI270731v1_random": "1651", "chr22_KI270732v1_random": "1652", "chr22_KI270733v1_random": "1653", "chr22_KI270734v1_random": "1654", "chr22_KI270735v1_random": "1655", "chr22_KI270736v1_random": "1656", "chr22_KI270737v1_random": "1657", "chr22_KI270738v1_random": "1658", "chr22_KI270739v1_random": "1659", "chr22_KI270875v1_alt": "1660", "chr22_KI270876v1_alt": "1661", "chr22_KI270877v1_alt": "1662", "chr22_KI270878v1_alt": "1663", "chr22_KI270879v1_alt": "1664", "chr22_KI270928v1_alt": "1665", "chr22_KN196485v1_alt": "1666", "chr22_KN196486v1_alt": "1667", "chr22_KQ458387v1_alt": "1668", "chr22_KQ458388v1_alt": "1669", "chr22_KQ759761v1_alt": "1670", "chr22_KQ759762v1_fix": "1671", "chr22_KQ759762v2_fix": "1672", "chr22_ML143378v1_fix": "1673", "chr22_ML143379v1_fix": "1674", "chr22_ML143380v1_fix": "1675", "chrUn_GL000195v1": "1676", "chrUn_GL000213v1": "1677", "chrUn_GL000214v1": "1678", "chrUn_GL000216v2": "1679", "chrUn_GL000218v1": "1680", "chrUn_GL000219v1": "1681", "chrUn_GL000220v1": "1682", "chrUn_GL000224v1": "1683", "chrUn_GL000226v1": "1684", "chrUn_KI270302v1": "1685", "chrUn_KI270303v1": "1686", "chrUn_KI270304v1": "1687", "chrUn_KI270305v1": "1688", "chrUn_KI270310v1": "1689", "chrUn_KI270311v1": "1690", "chrUn_KI270312v1": "1691", "chrUn_KI270315v1": "1692", "chrUn_KI270316v1": "1693", "chrUn_KI270317v1": "1694", "chrUn_KI270320v1": "1695", "chrUn_KI270322v1": "1696", "chrUn_KI270329v1": "1697", "chrUn_KI270330v1": "1698", "chrUn_KI270333v1": "1699", "chrUn_KI270334v1": "1700", "chrUn_KI270335v1": "1701", "chrUn_KI270336v1": "1702", "chrUn_KI270337v1": "1703", "chrUn_KI270338v1": "1704", "chrUn_KI270340v1": "1705", "chrUn_KI270362v1": "1706", "chrUn_KI270363v1": "1707", "chrUn_KI270364v1": "1708", "chrUn_KI270366v1": "1709", "chrUn_KI270371v1": "1710", "chrUn_KI270372v1": "1711", "chrUn_KI270373v1": "1712", "chrUn_KI270374v1": "1713", "chrUn_KI270375v1": "1714", "chrUn_KI270376v1": "1715", "chrUn_KI270378v1": "1716", "chrUn_KI270379v1": "1717", "chrUn_KI270381v1": "1718", "chrUn_KI270382v1": "1719", "chrUn_KI270383v1": "1720", "chrUn_KI270384v1": "1721", "chrUn_KI270385v1": "1722", "chrUn_KI270386v1": "1723", "chrUn_KI270387v1": "1724", "chrUn_KI270388v1": "1725", "chrUn_KI270389v1": "1726", "chrUn_KI270390v1": "1727", "chrUn_KI270391v1": "1728", "chrUn_KI270392v1": "1729", "chrUn_KI270393v1": "1730", "chrUn_KI270394v1": "1731", "chrUn_KI270395v1": "1732", "chrUn_KI270396v1": "1733", "chrUn_KI270411v1": "1734", "chrUn_KI270412v1": "1735", "chrUn_KI270414v1": "1736", "chrUn_KI270417v1": "1737", "chrUn_KI270418v1": "1738", "chrUn_KI270419v1": "1739", "chrUn_KI270420v1": "1740", "chrUn_KI270422v1": "1741", "chrUn_KI270423v1": "1742", "chrUn_KI270424v1": "1743", "chrUn_KI270425v1": "1744", "chrUn_KI270429v1": "1745", "chrUn_KI270435v1": "1746", "chrUn_KI270438v1": "1747", "chrUn_KI270442v1": "1748", "chrUn_KI270448v1": "1749", "chrUn_KI270465v1": "1750", "chrUn_KI270466v1": "1751", "chrUn_KI270467v1": "1752", "chrUn_KI270468v1": "1753", "chrUn_KI270507v1": "1754", "chrUn_KI270508v1": "1755", "chrUn_KI270509v1": "1756", "chrUn_KI270510v1": "1757", "chrUn_KI270511v1": "1758", "chrUn_KI270512v1": "1759", "chrUn_KI270515v1": "1760", "chrUn_KI270516v1": "1761", "chrUn_KI270517v1": "1762", "chrUn_KI270518v1": "1763", "chrUn_KI270519v1": "1764", "chrUn_KI270521v1": "1765", "chrUn_KI270522v1": "1766", "chrUn_KI270528v1": "1767", "chrUn_KI270529v1": "1768", "chrUn_KI270530v1": "1769", "chrUn_KI270538v1": "1770", "chrUn_KI270539v1": "1771", "chrUn_KI270544v1": "1772", "chrUn_KI270548v1": "1773", "chrUn_KI270579v1": "1774", "chrUn_KI270580v1": "1775", "chrUn_KI270581v1": "1776", "chrUn_KI270582v1": "1777", "chrUn_KI270583v1": "1778", "chrUn_KI270584v1": "1779", "chrUn_KI270587v1": "1780", "chrUn_KI270588v1": "1781", "chrUn_KI270589v1": "1782", "chrUn_KI270590v1": "1783", "chrUn_KI270591v1": "1784", "chrUn_KI270593v1": "1785", "chrUn_KI270741v1": "1786", "chrUn_KI270742v1": "1787", "chrUn_KI270743v1": "1788", "chrUn_KI270744v1": "1789", "chrUn_KI270745v1": "1790", "chrUn_KI270746v1": "1791", "chrUn_KI270747v1": "1792", "chrUn_KI270748v1": "1793", "chrUn_KI270749v1": "1794", "chrUn_KI270750v1": "1795", "chrUn_KI270751v1": "1796", "chrUn_KI270752v1": "1797", "chrUn_KI270753v1": "1798", "chrUn_KI270754v1": "1799", "chrUn_KI270755v1": "1800", "chrUn_KI270756v1": "1801", "chrUn_KI270757v1": "1802", "chrX_KI270880v1_alt": "1803", "chrX_KI270881v1_alt": "1804", "chrX_KI270913v1_alt": "1805", "chrX_KV766199v1_alt": "1806", "chrX_ML143381v1_fix": "1807", "chrX_ML143382v1_fix": "1808", "chrX_ML143383v1_fix": "1809", "chrX_ML143384v1_fix": "1810", "chrX_ML143385v1_fix": "1811", "chrX_MU273393v1_fix": "1812", "chrX_MU273394v1_fix": "1813", "chrX_MU273395v1_alt": "1814", "chrX_MU273396v1_alt": "1815", "chrX_MU273397v1_alt": "1816", "chrY_KI270740v1_random": "1817", "chrY_KN196487v1_fix": "1818", "chrY_KZ208923v1_fix": "1819", "chrY_KZ208924v1_fix": "1820", "chrY_MU273398v1_fix": "1821"}

    try:
        return ints[chrom]
    except Exception as e:
        pass # handel the exception later

def decipherVariant(chrom:str, pos:int, ref:str, alt:str) -> int:

    chrom = getChromINT(chrom=chrom)
    ref= getBaseINT(ref)
    alt = getBaseINT(alt)
    try:
        return int(chrom + str(pos) + ref + alt)
    except Exception as e:
        pass # handel the exception later


#######################

class FastaVCF:
    def __init__(self, fasta):
        self.fasta = fasta
        self.readFasta()


    def readFasta(self):

        Alternatives = {
            "A": ["C", "G", "T"],
            "C": ["A", "G", "T"],
            "G": ["A", "C", "T"],
            "T": ["A", "C", "G"]
        }
        NCtCommon = ["A", "C", "G", "T"]
        NCtids = {
            "R": ["G", "A"],
            "Y": ["C", "T"],
            "K": ["G", "T"],
            "M": ["A", "C"],
            "S": ["G", "C"],
            "W":["A", "T"],
        }
        output = None
        # output.write(self.vcfHeader)
        with open(self.fasta) as file:
            gfasta = file.read()
            
        location = 0    
        gfastasections = gfasta.strip().split(">")        

        for section in gfastasections[60:]:
            chrlines = section.split("\n")
            chro = chrlines[0]

            print(f"Processing {chro}")

            chrbases = chrlines[1:]
            chrbasesjoined = "".join(chrbases)
            chrbasesjoined = chrbasesjoined.upper()

            if chro:

                outputLocations = f"/home/genwork2/Mert/DecipherUCSCFasta/Deciphered_UCSC_REFGEN/{chro}_GenomicVariants.vcf"
                output = open(outputLocations, "w")
                output.write("##fileformat=VCFv4.2\n")
                output.write("##FORMAT=<ID=DP,Number=1,Type=Integer,Description='Approximate read depth (reads with MQ=255 or with bad mates are filtered)'>\n")
                output.write("##INFO=<ID=ID,Number=A,Type=Integer,Description='Allele count in genotypes, for each ALT allele, in the same order as listed'>\n")
                output.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tGENOMICVARIANTS\n")

                
                location = 0
                
                for base in chrbasesjoined:
                    location += 1
                    if base in NCtCommon:
                            
                        for Alt in Alternatives[base]:
                            try:
                                id = decipherVariant(chro, str(location), base, Alt)
                                if id != 1:
                                    output.write(f"{chro}\t{location}\t.\t{base}\t{Alt}\t1\t.\tID={id}\tDP\t1\n")
                                else:
                                    print("ERROR: One variant was missed!")
                                    exit()
                            except KeyError:
                                pass

                    elif base in NCtids:
                        for base in NCtids[base]:
                            for Alt in Alternatives[base]:
                                try:
                                    id = decipherVariant(chro, str(location), base, Alt)
                                    if id != 1:
                                        output.write(f"{chro}\t{location}\t.\t{base}\t{Alt}\t1\t.\tID={id}\tDP\t1\n")
                                    else:
                                        print("ERROR: One variant was missed!")
                                        exit()
                                except KeyError:
                                    pass

        if output:
            output.close()     


if __name__ == "__main__":
    fasta = "/home/genwork2/Mert/DecipherUCSCFasta/hg38.fa"
    t = FastaVCF(fasta)





