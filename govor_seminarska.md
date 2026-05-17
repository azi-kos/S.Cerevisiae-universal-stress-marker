# Govor — Univerzalni stresni markerji pri *S. cerevisiae*

*Trajanje: ~10 minut  |  Tempo: ~140 besed/min*

---

## Slajd 1 — Naslovna stran *(~25 s)*

Pozdravljeni. Danes vam bom predstavil seminarsko nalogo pri predmetu funkcijska genomika, v kateri sem se ukvarjal z analizo genske ekspresije pri kvasovki *Saccharomyces cerevisiae*. Naslov naloge je **Univerzalni stresni markerji** — cilj je bil najti gene, ki se pri kvasovki odzovejo na stres ne glede na to, kakšen tip stresa celico zadene.

---

## Slajd 2 — Uvod in cilji *(~55 s)*

Najprej nekaj besed o motivaciji. Analiza genske ekspresije nam omogoča, da pogledamo, kako se izražanje genov spreminja, ko se spremenijo okoljski pogoji. Vsaka celica ima približno šest tisoč genov, in vsak gen se lahko aktivira močneje ali šibkeje — to merimo z RNA-sekvenciranjem.

V svoji nalogi sem se osredotočil na **tri različne stresne pogoje**: hipotonično okolje, kjer celica zaradi premajhnega osmotskega tlaka nabrekne; stradanje aminokislin, kjer celici manjkajo gradniki za sintezo proteinov; in stradanje glukoze, kjer celica izgubi glavni vir energije. Za vsak pogoj sem imel tudi kontrolne celice v normalnih razmerah.

Glavno vprašanje raziskave je bilo: **kateri geni se diferencialno izražajo v vseh treh pogojih hkrati?** Če najdemo take gene, ti verjetno predstavljajo del splošnega stresnega odziva — torej univerzalni odgovor celice, neodvisen od konkretnega tipa stresorja.

---

## Slajd 3 — Podatki *(~55 s)*

Podatke sem pridobil iz javne baze GEO, konkretno iz dataseta **GSE201386**. Posebnost teh podatkov je, da gre za **single-cell RNA-sekvenciranje** — to pomeni, da je vsak vzorec ena posamezna celica, ne celotna kultura. Skupno sem imel 117 celic in v vsaki izmerjeno aktivnost približno 6000 genov.

Organizem je *Saccharomyces cerevisiae*, sev S288C — torej navadni pekovski kvas, ki je standardni modelni organizem v molekularni biologiji.

Razporeditev celic po pogojih je bila neenakomerna: največ celic, 47, je bilo iz pogoja stradanja glukoze, najmanj — le 12 — pa iz hipotoničnega okolja. To je pomembno omeniti, ker manjše število celic neposredno vpliva na statistično moč analize.

Iz teh 117 ločenih datotek sem nato sestavil eno samo ekspresijsko matriko velikosti približno 6000 genov krat 117 celic, ki je bila vhod za nadaljnjo analizo.

---

## Slajd 4 — Metode: DESeq2 *(~65 s)*

Za diferencialno ekspresijsko analizo sem uporabil orodje **DESeq2**, ki je standardno orodje na tem področju. Postopek je imel štiri korake.

Najprej sestava ekspresijske matrike. Nato sem za vsak stresni pogoj posebej s pomočjo DESeq2 statistično testiral, ali se izražanje vsakega gena pomembno razlikuje glede na kontrolo. Tretji korak je bilo filtriranje rezultatov — obdržal sem samo gene s prilagojeno p-vrednostjo manjšo od 0,05 in z absolutno vrednostjo log2FC večjo od ena, kar pomeni vsaj dvakratno spremembo v izražanju. In končno, četrti korak — preseki — torej iskanje genov, ki so signifikantni v vseh treh pogojih hkrati.

DESeq2 za vsak gen vrne dve glavni vrednosti: **log2FC**, ki pove velikost spremembe izražanja — log2FC enako 1 pomeni dvakratno povečanje, minus 1 pa polovično zmanjšanje izražanja. In **padj** — to je p-vrednost po Benjamini-Hochberg korekciji za večkratno testiranje, ki nam pove statistično zanesljivost rezultata.

---

## Slajd 5 — Rezultati po pogojih *(~55 s)*

Rezultati po posameznih pogojih so bili sledeči. Pri hipotoničnem stresu se je diferencialno izražalo **383 genov**, pri stradanju aminokislin **2487 genov**, pri stradanju glukoze pa kar **3056 genov**. Vidimo torej, da je odziv na stradanje veliko obsežnejši od osmotskega odziva.

Najpomembnejši rezultat — to, kar sem dejansko iskal — je presek. Skupno **150 genov** je bilo signifikantno diferencialno izraženih v vseh treh stresnih pogojih hkrati. Na desni strani vidite Vennov diagram, ki to lepo ilustrira: vsak krog predstavlja en pogoj, presek vseh treh pa je naših 150 univerzalnih stresnih markerjev.

To so kandidatni geni splošnega stresnega odziva pri kvasovki.

---

## Slajd 6 — Heatmap *(~45 s)*

Tukaj imamo prikaz teh 150 skupnih DEG-ov v obliki heatmapa — natančneje top 50 najmočneje reguliranih. Vrstice predstavljajo posamezne gene, stolpci tri stresne pogoje, barva pa vrednost log2FC — rdeča pomeni povečano izražanje, modra pa zmanjšano.

Ključno opažanje je, da je **smer regulacije za večino genov enaka v vseh treh pogojih**. To pomeni, da gen, ki je v hipotoničnem stresu povišan, je povišan tudi pri stradanju aminokislin in pri stradanju glukoze. Točno to bi pričakovali, če obstaja resnično usklajen splošen stresni odziv — in te podatki to potrjujejo.

---

## Slajd 7 — GO enrichment *(~60 s)*

Naslednji korak je bila funkcijska interpretacija. Z **GO enrichment analizo** sem preveril, ali so med naših 150 skupnih genov prekomerno zastopane kakšne funkcijske kategorije — torej, ali se ti geni statistično pogosteje pojavljajo v določenih bioloških procesih kot pri naključnem vzorcu.

Identificiral sem štiri obogatene procese. Prvi je **procesiranje mRNA** — celica torej selektivno regulira, katere RNA se predelujejo. Drugi je **DNA rekombinacija** — kar je pričakovano, saj stres povzroča poškodbe DNA. Tretji je **DNA biosintezni proces**, vezan na metabolizem nukleotidov med stresom. In četrti je **respiratorni metabolizem** — kar nakazuje preklop iz fermentacije v dihanje.

Pomembno opozorilo: ker je 150 genov razmeroma majhen vzorec, po BH korekciji rezultati niso strogo statistično značilni. Zato te GO rezultate obravnavam kot **eksploratorne** — kot smiselno biološko hipotezo, ne kot dokončno potrditev.

---

## Slajd 8 — Top up-regulated geni *(~75 s)*

Poglejmo zdaj posamezne gene z največjim povečanjem transkripcije.

Na prvem mestu je **ZNF1**, oziroma gen YFL052W. To je transkripcijski faktor, ki igra ključno vlogo pri presnovnem preklopu iz fermentacije v respiracijo. Ko je ZNF1 aktiviran, celica zmanjša prioriteto hitre rasti in biosinteze ribosomov, in se usmeri v energijsko bolj učinkovito stanje. Ker gre za transkripcijski faktor, eno povišanje ZNF1 povzroči kaskadne spremembe v stotih drugih genih — to lepo razloži, zakaj smo opazili tako veliko skupnih DEG-ov.

Drugi je **SME1**, gen YOR159C. To je sestavni del **spliceosoma** — kompleksa, ki izrezuje introne iz pre-mRNA. Pri kvasovkah introne sicer vsebuje le okoli 5 odstotkov genov, a so to ravno najmočneje izraženi geni, predvsem ribosomalni proteini. Povišana ekspresija SME1 torej odraža reorganizacijo splice-aparata, da lahko celica med stresom selektivno procesira ključne transkripte.

Tretji je **FMP25**, oziroma gen YLR077W. Ta protein sodeluje pri sestavljanju respiratornega kompleksa III v mitohondrijski elektronski transportni verigi. Med stresom celica potrebuje več mitohondrijske aktivnosti za proizvodnjo ATP, hkrati pa stres mitohondrije poškoduje preko reaktivnih kisikovih spojin — FMP25 torej pomaga pri njihovem vzdrževanju in obnavljanju.

---

## Slajd 9 — DUS3 *(~45 s)*

Med geni z najbolj konsistentno **znižano** ekspresijo izstopa **DUS3**, gen YLR401C. DUS3 sodeluje pri posttranskripcijski modifikaciji tRNA — natančneje pretvarja uridine v dihidrouridine v variabilni regiji tRNA. Te modifikacije so ključne za pravilno zvitje, fleksibilnost in stabilnost tRNA, kar neposredno vpliva na učinkovitost translacije.

Biološki pomen je jasen: med stresom celica zmanjša globalno raven sinteze proteinov. Manj modificiranih tRNA pomeni počasnejšo translacijo, kar je v skladu s tem, kar bi pričakovali — energetsko potratni procesi se zatrejo, viri se preusmerijo v preživetje. To se popolnoma sklada s prej omenjenim padcem ribosomske biogeneze, ki ga sproži ZNF1.

---

## Slajd 10 — GeneMANIA mreža *(~55 s)*

Da bi rezultate dodatno validiral, sem izvedel še analizo z orodjem **GeneMANIA**, ki razkrije mrežo fizičnih in genetskih interakcij med izbranimi geni. Rezultati so podprli ugotovitve GO analize.

Najprej se je pokazala velika **spliceosomska gruča** okoli SME1, s šestimi sorodnimi proteini — SMX2, SMX3, SMD1, SMD2, SMD3 in SMB1 — kar potrjuje pomen reorganizacije splice-aparata med stresom.

Druga zanimivost je povezava med FMP25 in proteinom SRM1, ki sodeluje pri transportu makromolekul med jedrom in citoplazmo — kar je pri stresu kritično.

In tretje, identificirana je bila povezava z **RTG1** — transkripcijskim faktorjem retrogradnega odziva. Retrogradni odziv je signalna pot iz mitohondrija v jedro, ki se aktivira ravno ob mitohondrijskem stresu. To je še en neodvisen dokaz, da je mitohondrijska obremenitev osrednji element opaženega odziva.

---

## Slajd 11 — Sinteza *(~55 s)*

Če rezultate povzamem, je skupni imenovalec analiziranih genov **energetski in respiratorni stres** ob prehodu na oksidativni metabolizem.

ZNF1 deluje kot **transkripcijsko stikalo**, ki preklopi celico iz fermentacije v respiracijo in zmanjša ribosomsko biogenezo. SME1 omogoča **reorganizacijo spliceosoma** za selektivno procesiranje pre-mRNA stresnih genov. FMP25 zagotavlja **mitohondrijsko zmogljivost** in obnovo ob ROS poškodbah. In DUS3 svoj prispevek doda preko **upočasnitve translacije**.

Vsi štirje geni torej delujejo usklajeno znotraj iste zgodbe: kadar mitohondriji ne morejo zagotoviti dovolj ATP — bodisi zaradi neposredne poškodbe respiratornega kompleksa bodisi zaradi spremembe vira ogljika — se sproži usklajen splošen stresni odziv pri kvasovki.

---

## Slajd 12 — Omejitve *(~35 s)*

Še na hitro o omejitvah analize. Prvič, DESeq2 normalizacija predpostavlja stabilno skupno količino RNA med vzorci, kar pri stradanju glukoze morda ne drži v celoti. Drugič, kot že rečeno, 150 genov je premalo za zanesljive GO rezultate po BH korekciji. In tretjič, single-cell RNA-seq podatki imajo veliko biološko variabilnost med posameznimi celicami, kar zmanjšuje statistično moč — še posebej pri stradanju glukoze, kjer je variabilnost največja.

---

## Slajd 13 — Zaključek *(~50 s)*

Za zaključek tri ključne ugotovitve.

**Prvič:** identificirali smo 150 genov skupnega stresnega odgovora — gene, ki se diferencialno izražajo v vseh treh stresnih pogojih hkrati.

**Drugič:** skupni vzorec je jasen — preklop iz fermentacije v respiracijo, reorganizacija RNA-aparata in upočasnitev translacije.

In **tretjič:** štirje ključni geni — ZNF1, SME1, FMP25 in DUS3 — predstavljajo kandidatne markerje splošnega stresnega odziva pri *S. cerevisiae*, vredne nadaljnje eksperimentalne validacije.

Hvala za pozornost. Vesel bom vaših vprašanj.
