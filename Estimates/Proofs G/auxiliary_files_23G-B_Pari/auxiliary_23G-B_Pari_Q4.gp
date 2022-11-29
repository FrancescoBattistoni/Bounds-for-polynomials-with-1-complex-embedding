/*
auxiliary_23G-B_Pari_Q4.gp

     BORDER
     t = 1: (1+(1/2)*y*z*a*b)*(1+y)*(1+y*z)*(1-y*z*a)*(1+z*a)*(1-z*a*b)*(1+a*b+(a*b)^2)*(1+a)*(1-a*b)*(1+b)
     INNER 

  R(y,z,a,b) = (1+(1/2)*y*z*a*b)*(1+y)*(1+y*z)*(1-y*z*a)*(1+z*a)*(1-z*a*b)*(1+a*b+(a*b)^2)*(1+a)*(1-a*b)*(1+b);
*/

{res1z=
y^226 - 1167838459/31440080*y^225 + 80239489830391/158458003200*y^224 - 
    167111710237102237/88736481792000*y^223 - 1118461324533023372459/44723186823168000*y^222 + 
    11450399779309236027011/39753943842816000*y^221 - 
    3905007591851884837533601/17173703740096512000*y^220 - 
    4930401842475156421270358873/412168889762316288000*y^219 + 
    179914890409729681028383437409/3297351118098530304000*y^218 + 
    54938651895176571587021043168731/237409280503094181888000*y^217 - 
    316764562221332605029826468274897/135662446001768103936000*y^216 - 
    8057006676007924014881056392081407/17093468196222781095936000*y^215 + 
    36558997248866515086919002761695883891/615364855064020119453696000*y^214 - 
    76039053981049754089649718111892083029/703274120073165850804224000*y^213 - 
    93242396047093070089082111356193650455433/88612539129218897201332224000*y^212 + 
    768731841697518353681277953663976186219937/202542946581071765031616512000*y^211 + 
    85266774131740801234991556101949802889385649/6380102817303760598495920128000*y^210 - 
    5525678410435203019084969981849108460367432271/68054430051240113050623148032000*y^209 - 
    18172291937073988840299928137600318578653635659/163330632122976271321495555276800*y^208 + 
    8565630653080834979561427637705959343785076749067/6533225284919050852859822211072000*y^207 + 
    11659290639573296871427230515534405077833288118247/52265802279352406822878577688576000*y^206 - 
    256805068162356334101377168108456391591642924575123/14933086365529259092251022196736000*y^205 + 
    9678294666560777106530978197094967288274726399454383/836252836469638509166057243017216000*y^204 
    + 84817562424641125638045359627789494746727121455167793/446001512783807204888563862942515200*y^\
    203 - 301337300808815936684137498204327464767604515076836137/1163482207262105751883210077241344\
    000*y^202 - 16177432407103365629000828348819471317561419915304979011/89200302556761440977712772\
    58850304000*y^201 + 175388058521281715697899608576349370527111298340523267/48544382343815750191\
    952529299865600*y^200 + 103748376409097771780327251569619311615327588671465790568237/6850583236\
    359278667088340934797033472000*y^199 - 36271980159949372269175890518282473450917855825963161475\
    4177/9134110981812371556117787913062711296000*y^198 - 
    12270941820283032374044219193329961928935672848067489262141551/10960933178174845867341345495675\
    2535552000*y^197 + 64737372247371704415843514643077469036767361859046813957132141/1753749308507\
    97533877461527930804056883200*y^196 + 518168423534904179574613624404846182131273412652173344606\
    4675351/7014997234031901355098461117232162275328000*y^195 - 
    28056609031130870169910443313558483309482350668077562365269729497/93533296453758684734646148229\
    76216367104000*y^194 - 98284750552033268703208753957406784264176279352391864714105090109/224479\
    91148902084336315075575142919281049600*y^193 + 974993876826232436290446362669290349514867941766\
    7046610874961155837/448959822978041686726301511502858385620992000*y^192 + 
    42099441465496407013479612870899873785439852480156406700001144654491/17958392919121667469052060\
    46011433542483968000*y^191 - 492072156752556113641998785065100290479448226864359775054735538066\
    7/34702208539365541002999150647563933188096000*y^190 - 
    3274808895356716883396259588342271900185069737093579636278481358280889/287334286705946679504832\
    96736182936679743488000*y^189 + 138254368709135823319421554223637249952262060047071504514994931\
    87788423/16419102097482667400276169563533106674139136000*y^188 + 
    33230837978135974496873321050685119426680166910549933871953840422498241/65676408389930669601104\
    678254132426696556544000*y^187 - 24032307169049250997771498413022756763337610081530016580273769\
    8420381509/52541126711944535680883742603305941357245235200*y^186 - 
    271284481941094756351207811675781364905130532051061331183423964401524521/1313528167798613392022\
    09356508264853393113088000*y^185 + 167848311841695453005130544593037124086069687701928723001158\
    293920827079453/7355757739672234995323723964462831790014332928000*y^184 + 
    114982071418376925420449735065774357010539386882542937591330246458328165821/1471151547934446999\
    0647447928925663580028665856000*y^183 - 3083574714671497523679929121605797564912157265941539330\
    523701342535650643219/29423030958688939981294895857851327160057331712000*y^182 - 
    813558318049401528573899933130228639347342590552754727355807802776450747481/2942303095868893998\
    1294895857851327160057331712000*y^181 + 3480588072550011875122996362259715497429705731073596613\
    85466749870407585391/784614158898371732834530556209368724268195512320*y^180 + 
    21584284803927433193215916270248573663289947340589192785555202512970274360739/23538424766951151\
    9850359166862810617280458653696000*y^179 - 2121065680632249322862544420781125303814169406738443\
    195315983720102804058539/1225959623278705832553953994077138631669055488000*y^178 - 
    88050342334013836865548024297048207449619121853275283733202356816848754543393/31384566355934869\
    3133812222483747489707278204928000*y^177 + 5844739630731565765012119697978945938458981797461634\
    864015040226673667356818563/941536990678046079401436667451242469121834614784000*y^176 + 
    99701208091443730465751208258890030225554041135087590172554143711893218442839/13450528438257801\
    1343062381064463209874547802112000*y^175 - 1538717231412251620765038239402411287328220515703912\
    0010837784470182629340847459/753229592542436863521149333960993975297467691827200*y^174 - 
    3912325942833299229904760414589463806702723656724920612031391990918499115313159/301291837016974\
    7454084597335843975901189870767308800*y^173 + 2641597315224237547972140310005303816555366585890\
    50224631397682230103932862525589/4304169100242496362977996194062822715985529667584000*y^172 - 
    52999409507273480818016068630680485640142591754226542056063736341943945585377219/30129183701697\
    474540845973358439759011898707673088000*y^171 - 50396702375495613663999304444276914028415724065\
    48679730836048639443015912005666721/30129183701697474540845973358439759011898707673088000*y^170 
    + 497596152256837109318656078721153743277837526466475880949072074493902846739559567/15064591850\
    848737270422986679219879505949353836544000*y^169 + 
    1373851341581539507867518658892796081317818347598967599881670101754963864801877843/334768707796\
    6386060093997039826639890210967519232000*y^168 - 
    160030387614425771852283640074994224334086684789669998573567246350270372476768167/7532295925424\
    36863521149333960993975297467691827200*y^167 - 589469649274602336197194733257505325535926695877\
    980462244664985337820390122602701/654982254384727707409695073009559978519537123328000*y^166 + 
    7691748425607634618989644041212005250502142752326654482811649029721023718967806229/753229592542\
    4368635211493339609939752974676918272000*y^165 + 
    4434089619160570678222650134505436705295136660003247212979589488974371086518746041/251076530847\
    4789545070497779869979917658225639424000*y^164 - 
    4407170896534689488403551225306700950664840668421558097520455030230438482416846849/107604227506\
    0624090744499048515705678996382416896000*y^163 - 
    771927793815863561356077753505790907612358394391931607319342205359577538197749019/2410334696135\
    79796326767786867518072095189661384704*y^162 + 213458240800328599277461260992566013693097430516\
    296889281351118669174700074328445671/15064591850848737270422986679219879505949353836544000*y^16\
    1 + 59762381802709459361208006230630499072816301718619197802304100932853561588826954881/1004306\
    1233899158180281991119479919670632902557696000*y^160 - 
    35718743017108062858777559650766830998502238892530765944749954173167007013662466993/83692176949\
    1596515023499259956659972552741879808000*y^159 - 
    97243865678815214614199744712583916603915906126415017793643773914801257499019426691/75322959254\
    24368635211493339609939752974676918272000*y^158 + 
    60021629197894766462950670859278238572741656241952693760890772029069713821272281543/53802113753\
    0312045372249524257852839498191208448000*y^157 + 
    92339777362127055931032461148067458733017718098007549989631839303710078044559699501/30129183701\
    69747454084597335843975901189870767308800*y^156 - 
    3759106108286593308964073767370882917528061257625095347126709239530428359121968567229/150645918\
    50848737270422986679219879505949353836544000*y^155 - 
    114552277465785755073760591651014401718703941442276221332233543682289791157342559079/2008612246\
    779831636056398223895983934126580511539200*y^154 + 
    3499979395063304129723864696484362045273508124387868870192153432897531139630744924329/753229592\
    5424368635211493339609939752974676918272000*y^153 - 
    4218205624985869679271830368496879655797146439161851153911236587890166719151777949/100430612338\
    99158180281991119479919670632902557696000*y^152 - 
    150823705219506012978898900128758666042003153912047295299024631723309462179312000669/2183274181\
    28242569136565024336519992839845707776000*y^151 + 
    774827042627994609653166954792473988634520012229324951466414436469248470188522896629/1255382654\
    237394772535248889934989958829112819712000*y^150 + 
    268535065324793308657677203263057737544661959819767563269958611517269953450043840711/3274911271\
    92363853704847536504779989259768561664000*y^149 - 
    6763695575225240238898207639724945447905641446521497245955309711853734967875480070813/200861224\
    6779831636056398223895983934126580511539200*y^148 - 
    7650250889020869634866979442193959009445628603114977003801099817119342835436887673017/602583674\
    0339494908169194671687951802379741534617600*y^147 + 
    90473873740542973520441694019928595398872127965510050817841339182516966532152362337303/75322959\
    25424368635211493339609939752974676918272000*y^146 + 
    75614533874256577700900243809434153577639795048245220783901283436963802722937386831999/15064591\
    850848737270422986679219879505949353836544000*y^145 - 
    18181222491554258204263481052358622079864351795453739310280584014998166602908760556303/55794784\
    6327731010015666173304439981701827919872000*y^144 - 
    162983720883737968449160601655197339429515540285271335947866941836425670219628286790031/7532295\
    925424368635211493339609939752974676918272000*y^143 + 
    102677606324901143357375716188428846324963403959106216298442758555550981212968739970333/1506459\
    185084873727042298667921987950594935383654400*y^142 + 
    17395548713323830815185917733885824422242119064558240856559204057180552116411744007757/25107653\
    0847478954507049777986997991765822563942400*y^141 - 
    382809552193443643417312497110591729622290961831505640675049189480261060750392579630749/3766147\
    962712184317605746669804969876487338459136000*y^140 - 
    295916279104014238025113250292100041530790700780713750323872561523581151737653650507213/1883073\
    981356092158802873334902484938243669229568000*y^139 + 
    1048833362779213227512414655173090416995818888558550088824543931772060597497261077541097/150645\
    91850848737270422986679219879505949353836544000*y^138 + 
    1582805902581877296980879870824638871026420675727834982066101796901529523174610840377731/753229\
    5925424368635211493339609939752974676918272000*y^137 + 
    38534501593568399785799519197245355383711132303329830100631955652678112113799807160819/30744065\
    0017321168784142585290201622570394976256000*y^136 + 
    370524998150169170061519795735182520053008326056152719802482958090135608456120741033253/3766147\
    962712184317605746669804969876487338459136000*y^135 - 
    7060248295779866042081798096174345583447338408011769482261736252243466068268196122667/167384353\
    89831930300469985199133199451054837596160*y^134 - 
    358450183957157726527805511162270972332622760860461214464922832344874996492322157946489/2353842\
    47669511519850359166862810617280458653696000*y^133 + 
    87396809697330513471637751352948611029123996171811546350918444277574429410950709619511/75322959\
    25424368635211493339609939752974676918272000*y^132 + 
    12571348330990449773924138845730015902116108846745563181288973420675621033431033069878223/25107\
    65308474789545070497779869979917658225639424000*y^131 + 
    54543741965809242799375627871325534364929323396409672633686881627359760022322656203153471/15064\
    591850848737270422986679219879505949353836544000*y^130 - 
    7640633874280324105420922781702770326635727299746603393887482328697639653346430382025989/753229\
    592542436863521149333960993975297467691827200*y^129 - 
    9019732213658323687873981197889040036912488363014370229903610085533011921030436712600067/602583\
    674033949490816919467168795180237974153461760*y^128 + 
    1272127423431614225347231963484843941309468319951476154591759974444106299121506453627073/109163\
    709064121284568282512168259996419922853888000*y^127 + 
    6070728780730964913930984631257083222982922006067557846189159629604595384321676999538271/163745\
    563596181926852423768252389994629884280832000*y^126 + 
    62737351423338880732789484818494308809725496908018326577640245544735253874190774465858001/15064\
    591850848737270422986679219879505949353836544000*y^125 - 
    1819368720455617723314627109880140171933299354334175678516434933828304599307029950124412267/301\
    29183701697474540845973358439759011898707673088000*y^124 - 
    7308868511327836630362959200328755603130949566114467458898496995429210368578775810439013/133907\
    483118655442403759881593065595608438700769280*y^123 + 
    59152788367711708155358273526772941266158901004762067587955822600743933460979586072054849/13099\
    64508769455414819390146019119957039074246656000*y^122 + 
    284001980593962141636252004829255834909864256068234360702781075020009066220475529315028259/2152\
    084550121248181488998097031411357992764833792000*y^121 + 
    837938268738576387334632150564507631109772971292829208981425936606020862173950158395184981/1004\
    3061233899158180281991119479919670632902557696000*y^120 - 
    180159155078314609482267433351420211845851545666120436262587880763168178188685083669802409/1255\
    382654237394772535248889934989958829112819712000*y^119 - 
    91381700434341299153684234830737630026008426290092649548344072987399215655220506777564433/23912\
    0505569027575720999788559045706443640537088000*y^118 - 
    3215350695549745081751095269069668225694664423345050307514956100517409673670911298748921/219600\
    46429808654913153041806442973040742498304000*y^117 + 
    788678628826348721847850671742868439579460468995068391696524176409376266207867321790620963/1076\
    042275060624090744499048515705678996382416896000*y^116 + 
    7971845764465661802332373611527352362941338313601238157013080557460616417596812284856009563/753\
    2295925424368635211493339609939752974676918272000*y^115 - 
    18826296859012481145900445450474365855313006892044056290397035677627747784163195556153755851/30\
    129183701697474540845973358439759011898707673088000*y^114 - 
    40234159877941808819954602639418673574857292538702088072302754267022859859156189078340011293/15\
    064591850848737270422986679219879505949353836544000*y^113 - 
    5305561304644318094574374319656108867708577165787326109907453729069324381241730224189638077/602\
    5836740339494908169194671687951802379741534617600*y^112 + 
    31867518778109425183504030211189062624121842177019012119566971693091989365736619638367994691/75\
    32295925424368635211493339609939752974676918272000*y^111 + 
    1645563978458747076104834423704219347319089286020330644012924642504264844890122050536915359/358\
    680758353541363581499682838568559665460805632000*y^110 - 
    2845897349850581871994262769837700860470811078063146477003062976515497700604223280114411577/753\
    229592542436863521149333960993975297467691827200*y^109 - 
    148524947356764293488543676832131931568027363753155149263977925700031554202159539278162241083/1\
    5064591850848737270422986679219879505949353836544000*y^108 - 
    15798779873397522020931312036736049603821231888064674216040703792049843367632026617400407333/15\
    064591850848737270422986679219879505949353836544000*y^107 + 
    404270170453217079785186169547074790509874156079118433750806415019037044291753231639789355811/3\
    0129183701697474540845973358439759011898707673088000*y^106 + 
    78610643339509899512315040312313599584643136949450726905966816634399804908208859158218088681/75\
    32295925424368635211493339609939752974676918272000*y^105 - 
    62465056703106664422592920107161539819527214633699683215090656318819478332613572536193200639/60\
    25836740339494908169194671687951802379741534617600*y^104 - 
    4271086480653546356962102053744463114400858284999255714615814441917733312228275435007605597/218\
    327418128242569136565024336519992839845707776000*y^103 - 
    461278616132549928241961503618811143273983010445866222512008101796765812699515930061936201/4707\
    68495339023039700718333725621234560917307392000*y^102 + 
    18723056712956036173827062515952441771534921310786086519149328530448638773292403712395822737/94\
    1536990678046079401436667451242469121834614784000*y^101 + 
    79149267790798652007593229216241415899049375886924541001693485436235380136142833249262553431/60\
    25836740339494908169194671687951802379741534617600*y^100 - 
    251249596647753464545740921597902311313261713285136123458863953576629191123068743848445520107/3\
    0129183701697474540845973358439759011898707673088000*y^99 - 
    1060982400011692329504467175613181281908530148379086003965981818895023630707274984578980813/935\
    68893483532529629956439001365711217076731904000*y^98 - 
    288320737704696702497977792448820970719243301812116684122326581729542341020372793360142411/6025\
    83674033949490816919467168795180237974153461760*y^97 - 
    10450866623895800030621684423772129168498834972160013751325302153161773869434476462019103647/12\
    55382654237394772535248889934989958829112819712000*y^96 - 
    24387668269132049085187797550090558286994880437411248223880694598045274847271015595908799847/10\
    76042275060624090744499048515705678996382416896000*y^95 + 
    51490470018264396683070844792817178606601370899063211561846017742192242172754765918144725611/37\
    66147962712184317605746669804969876487338459136000*y^94 + 
    161108406064321104585009011423721618217567549142974480100680137876566354669399604462135666287/1\
    883073981356092158802873334902484938243669229568000*y^93 + 
    12466770191756132133600414919658427733204663608267423488700271486801891787367967609260718591/20\
    9230442372899128755874814989164993138185469952000*y^92 - 
    45997666254194241983338043732913177671576180115965975540310392012450856200341502298351920627/37\
    6614796271218431760574666980496987648733845913600*y^91 - 
    24904572178220342927964736370368728426358319358628361528860409707095886745661359433063563747/10\
    4615221186449564377937407494582496569092734976000*y^90 - 
    5137962634211769482939708712758752058643605143755472622374634124657962654789634379708727029/941\
    536990678046079401436667451242469121834614784000*y^89 + 
    2075065038237862742873491825956882459738849919253979530909587282152735261004963242798832627/511\
    7048862380685214138242757887187332183883776000*y^88 + 
    59508592186187445434955674440562859113896067682800996327451172976357650131536745725417585657/15\
    6922831779674346566906111241873744853639102464000*y^87 - 
    69431408897838024969200931956574483275904060045452261286434161030153877075251222297548392539/23\
    5384247669511519850359166862810617280458653696000*y^86 - 
    194668361968391367554441684551206656118941067928828903105945418794079980993818950025713728411/2\
    35384247669511519850359166862810617280458653696000*y^85 - 
    3822253353587543424701795371650973091316143820880358610103991440279478261886125356102458689/130\
    76902648306195547242175936822812071136591872000*y^84 + 
    10491645437033094524171992787142855360091171974725278127969547745490101547038213072324534313/11\
    769212383475575992517958343140530864022932684800*y^83 + 
    9593736042525747528494927663551825883340376613015522464393841159815243941937966317918657841/840\
    6580273911125708941398816528950617159237632000*y^82 - 
    2161023043397123427916682212061675051350882847981458974870723919853964244369570270282243547/117\
    69212383475575992517958343140530864022932684800*y^81 - 
    45799321018505827175507915931871713049277179639353343659757087582608816860771337208663881271/29\
    423030958688939981294895857851327160057331712000*y^80 - 
    108877269323250352570134069239595026500350953013422542718901378553531637461699849505392513/1021\
    63301939892152712829499506428219305754624000*y^79 + 
    13307425340973986224455431801613710513832113721368529008852589950520444076838390470302094613/14\
    711515479344469990647447928925663580028665856000*y^78 + 
    1756189593731090313774550954972019485518982296390568769283553264952795751519955900412345803/919\
    469717459029374415465495557853973751791616000*y^77 + 
    2168573093023784850914716831011802881796145833694121595153966509941597423650625742466979707/367\
    7878869836117497661861982231415895007166464000*y^76 - 
    1355662015510866788685391594582470609817062642287004942614383868549508999811194802127203553/919\
    469717459029374415465495557853973751791616000*y^75 - 
    648360263455317263471680821868185767486202018259598917015882348673054888016635614021520407/3677\
    87886983611749766186198223141589500716646400*y^74 - 
    266652149798630430081064062627903300518895592139017596190912882466525183021270961388523/1998847\
    2118674551617727510772996825516343296000*y^73 + 14631616195185968789601053268212469637569165824\
    91289606868316133822391954604676259344707493/919469717459029374415465495557853973751791616000*y\
    ^72 + 908811344072740021099515121396554932377787274189989537214128761325479230321346322798423/7\
    18335716764866698762082418404573416993587200*y^71 - 
    79065055933576048838675904949467368553923343351880030065102813064273363828149924848227679/22986\
    7429364757343603866373889463493437947904000*y^70 - 
    779245147039495469324127623845659973830964509386461185441604865653260129326486490352741/6081148\
    92499358051862080354204929876819968000*y^69 - 5143023740540981847651635947472408171885330253713\
    505709421678676399456389999222724369043/7183357167648666987620824184045734169935872000*y^68 + 
    3843919641931940834791890591988702271810465201637593317880248357117315063287616643521229/957780\
    9556864889316827765578727645559914496000*y^67 + 25033357917475229428044472951574212097799383210\
    6198041952613642843079103500032319464521/312319876854289869026992355828075398692864000*y^66 + 
    227668983172351289309919474694826101988028144886546007883493698854995176087391219980107/7183357\
    16764866698762082418404573416993587200*y^65 - 3542697639141233948119836691218301270387084310181\
    90872908764699555729213321302592223/1239364590691626464392826808841569042432000*y^64 - 
    118149642002932359819177923132432602532142861523150173960592179073786057967964638751749/2993065\
    48652027791150867674335238923747328000*y^63 - 9439817580962293384941734811294459695534994903613\
    8047173870257820626869776939448536301/897919645956083373452603023005716771241984000*y^62 + 
    33308837876823368841610756437005062526476053284576144910742389051544492531101990663617/22447991\
    1489020843363150755751429192810496000*y^61 + 17123611788437714422531893531082822721084588305033\
    865238702747605983099500934306607863/112239955744510421681575377875714596405248000*y^60 + 
    837529314077508438199179190979973809211826908902139692092086309877802772167020704361/3741331858\
    1503473893858459291904865468416000*y^59 - 16597030084063081098070100467210595397805932132173822\
    40916168667559551272609443789929/28059988936127605420393844468928649101312000*y^58 - 
    23597376391056788565347798215149328014788849395069512654345045363915117094649009827/51962942474\
    3103804081367490165345353728000*y^57 - 90233215565276259499895263363924093271115363407374038983\
    6001242994367870071806441/2004284924009114672885274604923474935808000*y^56 + 
    1162175539007592593449641470056466932726290162887365541737269507000082630958440399/626339038752\
    84833527664831403858591744000*y^55 + 3538852635770600457690863184300867951048686974735224578498\
    381509526825287859616907/350749861701595067754923055861608113766400*y^54 - 
    3262752597700431323552414049192179267057387576583697147705662889042414091947933393/175374930850\
    7975338774615279308040568832000*y^53 - 20296411187591343231780339883837122364317572524697151411\
    09642992220148366277793701/438437327126993834693653819827010142208000*y^52 - 
    11658479415345442584290311104148736932604770750582280500780088685295658079968533/73072887854498\
    97244894230330450169036800*y^51 + 1421609498418344063026100307254598091962742973084609373234354\
    90590028005876113/173983066320235648687957865010718310400*y^50 + 
    12377798952742096427835203827349547975310797077361432758617502458655948936978447/13701166472718\
    557334176681869594066944000*y^49 + 544058133766505376736163852839687193615654700192249883309845\
    21733545906486563/342529161817963933354417046739851673600*y^48 - 
    673243162552725396115527428387757950687675080704355152437311591775564306463011/3425291618179639\
    333544170467398516736000*y^47 - 123734225446661949176266964581788166421556042179505893244987090\
    22541630429119/95146989393878870376226957427736576000*y^46 - 
    92774328374501293205712100441117691879789313021982496677001550166504077319/30582960876603922620\
    930093458915328000*y^45 + 329857259214331582294715209586660490186492122894300696587320938738644\
    3001311/107040363068113729173255327106203648000*y^44 + 
    14482900253998301815172188879464724457856887647256766686474562577434700643/11634822072621057518\
    83210077241344000*y^43 - 3156114048439306244079679650936706848678233770034941810621131703635792\
    71/136531075341981797414866488655872000*y^42 - 115935880550536434338205073872376998115769090944\
    60917839443684847938610701/3345011345878554036664228972068864000*y^41 - 
    296530093099693836625761992799773101547717283288083990023333477285425977/4181264182348192545830\
    28621508608000*y^40 + 9129766475072958502227622050228659769252665546361524740130218569188531/19\
    910781820705678789668029595648000*y^39 + 642648426739739948191745083259634198520716815535618374\
    4974525461111201/20906320911740962729151431075430400*y^38 + 
    161209670663176685376275338673797501628045492078414420900071888597517/4544852372117600593293789\
    364224000*y^37 - 891681676397831682557202047012128187228999887636118672953891409088679/26132901\
    139676203411439288844288000*y^36 - 158735502481545236031594224164825499815147321022316452250208\
    20680093/933317897845578693265688887296000*y^35 - 
    936930549620009073968870711955217738221478717326290271682457952857/3629569602732806029366567895\
    04000*y^34 - 164564094574114739570951362425817223130016225295313959264503010641/816653160614881\
    356607477776384000*y^33 - 776121043064264824610466181243077407580132793048981246932677089/12602\
    67223171113204641169408000*y^32 - 9323732287789667650215483710107623977317190124354048287106624\
    3923/204163290153720339151869444096000*y^31 - 3303687365315589945626186264320570490283808967277\
    805921124460083/51040822538430084787967361024000*y^30 + 
    765021087413116008546610131485681947440010035716427141216577279/6380102817303760598495920128000\
    *y^29 + 102973714035140612419097460716995626166174651865534512333176179/91144325961482294264227\
    4304000*y^28 + 1011253533166055365831282413279434336015115821526113772075757/196916753620486438\
    22518272000*y^27 + 2034027682512372489496296422114643553251074256175839468934499/19937821304074\
    2518702997504000*y^26 - 548551872962520514113635593350613031621122265052491208495209/1993782130\
    40742518702997504000*y^25 - 148149866009732957194285588645856910029134846670164893622811/498445\
    53260185629675749376000*y^24 - 139512754746000237495225173949030266337575513733869326121/103842\
    819292053395157811200*y^23 - 3687158598544672534702993508110712505330835826034373531159/6230569\
    157523203709468672000*y^22 - 503095787774031579882959940381311627205495390735027827999/15576422\
    89380800927367168000*y^21 - 27327090543634669205644383293293375556366166502149403357/1947052861\
    72600115920896000*y^20 - 188118923239674694022621383275628295704515487012249423/108169603429222\
    28662272000*y^19 + 13413844355102766114557320306540306848411575452250021/5290904515559785758720\
    00*y^18 + 64097538459125256744173951123248457697428556496838953/3042270096446876811264000*y^17 +
    688529944057310572969880602110799461747554029993183/84507502679079911424000*y^16 + 
    77799127664878045200877394859692475473395742588303/63380627009309933568000*y^15 - 
    14983994516237479716125171439432484765317140247841/31690313504654966784000*y^14 - 
    8855532553441639111552252335418406675107845884861/23767735128491225088000*y^13 - 
    39586501911690700741278210543058979658796366601/330107432340155904000*y^12 - 
    5477186179415905960172819076217024537897595711/297096689106140313600*y^11 + 
    15215374201309249576100786713354680225569099/13263245049381264000*y^10 + 
    429666021604879833865730036005477895921101/331581126234531600*y^9 + 
    2065178403680381211639960195969229360388/6907940129886075*y^8 + 
    867231321440910983960039846684294854816/31536248419045125*y^7 - 
    1195075733416169654523606364477278106624/725333713638037875*y^6 - 
    92536980596736084591744709214474985472/145066742727607575*y^5 - 
    25598726894225662975675427898020200448/725333713638037875*y^4 + 
    4997497066955396563516479508578304/1074568464648945*y^3 + 
    44976547482541855529141635920166912/80592634848670875*y^2 - 
    2790672082866463322969207239344128/725333713638037875*y - 
    429705988522244675928385216577536/241777904546012625
;}


{res1a=
y^13*z^10 - 2/5*y^13*z^8 + 23/20*y^12*z^10 + 183/40*y^12*z^9 - 37/20*y^12*z^7 - 37/16*y^11*z^10 + 
    97/20*y^11*z^9 + 751/80*y^11*z^8 + 11/80*y^11*z^7 - 581/80*y^11*z^6 - 343/128*y^10*z^13 + 
    343/320*y^10*z^11 - 1193/320*y^10*z^10 - 6673/640*y^10*z^9 + 639/80*y^10*z^8 + 4213/640*y^10*z^7
    + 247/160*y^10*z^6 - 8253/320*y^10*z^5 - 12397/1280*y^9*z^13 - 6321/1280*y^9*z^12 + 
    3381/1280*y^9*z^11 + 699/1280*y^9*z^10 - 19959/1280*y^9*z^9 - 23823/1280*y^9*z^8 + 27/20*y^9*z^7
    - 8321/640*y^9*z^6 + 10851/1280*y^9*z^5 - 68337/1280*y^9*z^4 - 4319/320*y^8*z^13 - 
    1057/64*y^8*z^12 - 623/1280*y^8*z^11 + 1343/640*y^8*z^10 + 309/64*y^8*z^9 - 15011/640*y^8*z^8 - 
    4009/640*y^8*z^7 - 2269/80*y^8*z^6 - 19821/1280*y^8*z^5 + 2997/128*y^8*z^4 - 9687/160*y^8*z^3 - 
    733/80*y^7*z^13 - 2669/128*y^7*z^12 - 8197/1280*y^7*z^11 + 17167/1280*y^7*z^10 + 
    16747/1280*y^7*z^9 + 5739/256*y^7*z^8 - 6293/1280*y^7*z^7 + 58171/1280*y^7*z^6 - 
    95151/1280*y^7*z^5 + 42179/1280*y^7*z^4 + 2175/64*y^7*z^3 - 151/4*y^7*z^2 - 61/20*y^6*z^13 - 
    201/16*y^6*z^12 - 4203/640*y^6*z^11 + 59553/1280*y^6*z^10 + 40469/1280*y^6*z^9 + 
    28713/640*y^6*z^8 + 3309/160*y^6*z^7 + 28199/640*y^6*z^6 + 117867/1280*y^6*z^5 - 
    133867/1280*y^6*z^4 + 4563/64*y^6*z^3 + 1053/40*y^6*z^2 - 61/5*y^6*z - 2/5*y^5*z^13 - 
    29/8*y^5*z^12 - 833/320*y^5*z^11 + 3717/64*y^5*z^10 + 31823/320*y^5*z^9 + 719/40*y^5*z^8 + 
    45221/1280*y^5*z^7 - 30059/1280*y^5*z^6 + 102297/1280*y^5*z^5 + 91033/1280*y^5*z^4 - 
    491/5*y^5*z^3 + 4033/80*y^5*z^2 + 103/10*y^5*z - 8/5*y^5 - 2/5*y^4*z^12 - 29/80*y^4*z^11 + 
    5299/160*y^4*z^10 + 72231/640*y^4*z^9 + 32449/1280*y^4*z^8 - 3707/256*y^4*z^7 - 
    14453/640*y^4*z^6 - 77517/1280*y^4*z^5 + 82667/1280*y^4*z^4 + 8107/320*y^4*z^3 - 
    9569/160*y^4*z^2 + 61/4*y^4*z + 8/5*y^4 + 1417/160*y^3*z^10 + 8767/160*y^3*z^9 + 
    3027/128*y^3*z^8 - 8687/160*y^3*z^7 - 3283/320*y^3*z^6 - 30459/640*y^3*z^5 - 15103/320*y^3*z^4 +
    501/16*y^3*z^3 + 1401/160*y^3*z^2 - 811/40*y^3*z + 8/5*y^3 + 9/10*y^2*z^10 + 1851/160*y^2*z^9 + 
    741/320*y^2*z^8 - 14841/320*y^2*z^7 - 36063/1280*y^2*z^6 + 807/80*y^2*z^5 - 5923/256*y^2*z^4 - 
    5333/320*y^2*z^3 + 1667/160*y^2*z^2 + 26/5*y^2*z - 14/5*y^2 + 33/40*y*z^9 - 459/160*y*z^8 - 
    3087/160*y*z^7 - 6999/320*y*z^6 + 2661/640*y*z^5 + 6861/640*y*z^4 - 205/64*y*z^3 - 529/160*y*z^2
    + 29/20*y*z + 13/10*y - 3/5*z^8 - 243/80*z^7 - 393/80*z^6 - 231/160*z^5 + 261/80*z^4 + 
    417/160*z^3 - 23/40*z - 1/5
;}


{res1b=
y^4*z^2*a - 5/2*y^3*z^3*a^2 + 3/4*y^3*z^2*a + 3/4*y^3*z*a - 3/4*y^3*z - 19/4*y^2*z^3*a^2 - 
    7/4*y^2*z^2*a^2 + 1/4*y^2*z^2*a + 1/2*y^2*z*a - 1/2*y^2*z - 1/2*y^2 - 2*y*z^3*a^2 - 3*y*z^2*a^2 
    + 2*y*z^2*a + y*z - 1/4*y - z^2*a^2 + z^2*a + 3/4*z*a + 1/2*z + 1/2
;}


{Ly=
y^3*z^3*a^2*b + 3/4*y^2*z^3*a^2*b + 3/4*y^2*z^2*a^2*b - 3/4*y^2*z^2*a*b + 3/2*y^2*z^2*a + 
    1/2*y*z^2*a^2*b - 1/2*y*z^2*a*b + y*z^2*a - 1/2*y*z*a*b + y*z*a - y*z - 1/4*z*a*b + 1/2*z*a - 
    1/2*z - 1/2
;}


