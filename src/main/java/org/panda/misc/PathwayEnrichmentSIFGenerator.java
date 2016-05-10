package org.panda.misc;

import org.biopax.paxtools.controller.Cloner;
import org.biopax.paxtools.controller.Completer;
import org.biopax.paxtools.controller.SimpleEditorMap;
import org.biopax.paxtools.io.SimpleIOHandler;
import org.biopax.paxtools.model.BioPAXElement;
import org.biopax.paxtools.model.BioPAXLevel;
import org.biopax.paxtools.model.Model;
import org.biopax.paxtools.model.level3.Pathway;
import org.biopax.paxtools.pattern.miner.SIFEnum;
import org.biopax.paxtools.pattern.miner.SIFInteraction;
import org.biopax.paxtools.pattern.miner.SIFSearcher;
import org.biopax.paxtools.pattern.miner.SIFType;
import org.biopax.paxtools.pattern.util.Blacklist;
import org.panda.utility.Kronometre;
import org.panda.utility.graph.SIFGenerator;

import java.awt.*;
import java.io.FileInputStream;
import java.io.IOException;
import java.util.*;
import java.util.function.Function;

import static java.util.stream.Collectors.toMap;
import static java.util.stream.Collectors.toSet;

/**
 * @author Ozgun Babur
 */
public class PathwayEnrichmentSIFGenerator extends SIFGenerator
{
	Set<String> genes;

	Set<String> pathwayIDs;

	String[][] groupTerms;

	String owlFilename;
	String blacklistFilename;

	SIFType[] types;

	private static final Color[] COLORS = new Color[]{
		new Color(0, 69, 134),
		new Color(255, 66, 14),
		new Color(255, 211, 32),
		new Color(87, 157, 28),
		new Color(126, 0, 33),
		new Color(131, 202, 255),
		new Color(49, 64, 4),
		new Color(174, 207, 0),
		new Color(75, 31, 111),
		new Color(255, 149, 14),
		new Color(197, 0, 11),
		new Color(0, 132, 209),
		new Color(255, 128, 128),
		new Color(119, 33, 111),
		new Color(102, 153, 153),
	};

	private Set<SIFInteraction> getSIF(Set<Pathway> pathways, Model model, Blacklist blacklist, SIFType... types)
	{
		model = excise(model, pathways);
		SIFSearcher searcher = new SIFSearcher(types);
		searcher.setBlacklist(blacklist);
		return searcher.searchSIF(model);
	}

	private static Model excise(Model model, Set<Pathway> pathways)
	{
		Completer c = new Completer(SimpleEditorMap.L3);
		Set<BioPAXElement> objects = c.complete(new ArrayList<>(pathways), model);
		Cloner cln = new Cloner(SimpleEditorMap.L3, BioPAXLevel.L3.getDefaultFactory());
		return cln.clone(model, objects);
	}

	private Set<SIFInteraction> removeIrrelevant(Set<SIFInteraction> sifs)
	{
		return sifs.stream().filter(i -> genes.contains(i.sourceID) || genes.contains(i.targetID))
			.collect(toSet());
	}

	protected void prepare() { try
	{
		setDefaultNodeBorderWidth(2);
		setDefaultNodeBorderColor(Color.BLACK);

		// read model and blacklist
		SimpleIOHandler h = new SimpleIOHandler(BioPAXLevel.L3);
		Model model = h.convertFromOWL(new FileInputStream(owlFilename));
		Blacklist blacklist = new Blacklist(blacklistFilename);

		// find pathways
		Set<Pathway> pathways = pathwayIDs.stream().map(id -> (Pathway) model.getByID(id)).collect(toSet());

		// groups pathways
		Set<Set<Pathway>> groups = Arrays.stream(groupTerms).map(terms -> select(pathways, terms))
			.collect(toSet());

		// get SIFs of pathway groups
		Map<Set<Pathway>, Set<SIFInteraction>> p2sif = groups.parallelStream().collect(toMap(
			Function.identity(), g -> removeIrrelevant(getSIF(g, model, blacklist, types))));

		// assign colors to pathway sets
		Map<Set<Pathway>, Color> colorMap = new HashMap<>();
		int i = 0;
		for (Set<Pathway> ps : p2sif.keySet())
		{
			colorMap.put(ps, COLORS[i++]);
		}
		printPathwayColors(colorMap);

		Set<String> linkerGenes = getGenesInMultiplePathwayGroups(p2sif.values());

		p2sif.forEach((ps, sifs) ->
		{
			sifs.forEach(s -> addEdge(s.sourceID, s.targetID, s.type.getTag(), s.getMediatorsInString()));
			Color color = colorMap.get(ps);
			extractGenes(sifs).stream().filter(gene -> !linkerGenes.contains(gene)).forEach(gene ->
				addNodeBorderColor(gene, color));
		});

		Set<String> genesInSIFs = p2sif.values().stream().flatMap(Collection::stream)
			.map(s -> new String[]{s.sourceID, s.targetID}).flatMap(Arrays::stream).collect(toSet());

		genesInSIFs.stream().filter(genes::contains).filter(g -> !hasNodeColor(g)).forEach(gene ->
			addNodeColor(gene, Color.LIGHT_GRAY));
	}
	catch (IOException e){throw new RuntimeException(e);}}

	/**
	 * This method is useful for creating a legend for the generated SIF graph. It prints the included pathway names in
	 * html style, which can be displayed in a markdown-supporting wiki page.
	 * @param colorMap
	 */
	private void printPathwayColors(Map<Set<Pathway>, Color> colorMap)
	{
		System.out.println("<html><body><b>");
		colorMap.keySet().forEach(set ->
		{
			Color c = colorMap.get(set);
			set.forEach(p -> System.out.println("<span style=\"color: rgb(" +
				c.getRed() + "," + c.getGreen() + "," + c.getBlue() + ")\">" + p.getDisplayName() + "</span><br>"));
		});
		System.out.println("</b></body></html>");
	}

	private Set<String> getGenesInMultiplePathwayGroups(Collection<Set<SIFInteraction>> sifGroups)
	{
		Map<String, Integer> cnt = new HashMap<>();

		sifGroups.forEach(sifs -> extractGenes(sifs).forEach(gene ->
				cnt.put(gene, cnt.containsKey(gene) ? cnt.get(gene) + 1 : 1)));

		return cnt.keySet().stream().filter(gene -> cnt.get(gene) > 1).collect(toSet());
	}

	private Set<String> extractGenes(Set<SIFInteraction> sifs)
	{
		return sifs.stream().map(s -> new String[]{s.sourceID, s.targetID}).flatMap(Arrays::stream).collect(toSet());
	}

	/**
	 * Gets the pathways that contains any of the given terms in its name.
	 */
	private Set<Pathway> select(Set<Pathway> pathways, String[] terms)
	{
		return pathways.stream().filter(p -> Arrays.stream(terms).anyMatch(p.getDisplayName().toLowerCase()::contains))
			.collect(toSet());
	}

	public void setGenes(Set<String> genes)
	{
		this.genes = genes;
	}

	public void setPathwayIDs(Set<String> pathwayIDs)
	{
		this.pathwayIDs = pathwayIDs;
	}

	public void setGroupTerms(String[][] groupTerms)
	{
		for (int i = 0; i < groupTerms.length; i++)
		{
			for (int j = 0; j < groupTerms[i].length; j++)
			{
				groupTerms[i][j] = groupTerms[i][j].toLowerCase();
			}
		}
		this.groupTerms = groupTerms;
	}

	public void setOwlFilename(String owlFilename)
	{
		this.owlFilename = owlFilename;
	}

	public void setBlacklistFilename(String blacklistFilename)
	{
		this.blacklistFilename = blacklistFilename;
	}

	public void setTypes(SIFType... types)
	{
		this.types = types;
	}

	public static void main(String[] args)
	{
		Kronometre k = new Kronometre();
		PathwayEnrichmentSIFGenerator gen = new PathwayEnrichmentSIFGenerator();
		gen.setOwlFilename("/home/babur/Documents/PC/PathwayCommons.8.Detailed.BIOPAX.owl");
		gen.setBlacklistFilename("/home/babur/Documents/PC/blacklist_v8.txt");
		Set<String> pathwayIDs = new HashSet<>(Arrays.asList(
			"http://pathwaycommons.org/pc2/Pathway_10f0f96341b68cece9ba6646e2222c56",
			"http://pathwaycommons.org/pc2/Pathway_1d37e4917ffea3abd6fefbe484d32787",
			"http://pathwaycommons.org/pc2/Pathway_7869cf22fdff2057312c151361d24bb7",
			"http://identifiers.org/panther.pathway/P00034",
			"http://pathwaycommons.org/pc2/Pathway_a1fad3aea6bed5c73670bb81501e3294",
			"http://identifiers.org/reactome/R-HSA-1650814",
			"http://identifiers.org/reactome/R-HSA-1474290"));
		gen.setPathwayIDs(pathwayIDs);
		gen.setGroupTerms(new String[][]{new String[]{"collagen"}, new String[]{"integrin"}});
		gen.setTypes(SIFEnum.CONTROLS_STATE_CHANGE_OF, SIFEnum.CONTROLS_EXPRESSION_OF);
		Set<String> genes = new HashSet<>(Arrays.asList("CDA, SERPINE2, HOXA9, GJA1, MPRIP, MYC, DPYSL3, PRSS3, PRSS2, IER3, GSTK1, HPCAL1, DKK1, FAR2, HPRT1, PRR13, GNE, HS3ST3A1, COL13A1, CRABP2, CDCA4, S100A16, RDH10, S100A13, SFN, PLEK2, S100A10, APCDD1L, NR2F1, FZD8, MLLT11, FOSL1, COL1A2, PKIA, TRIP6, FXYD5, PAM, FERMT2, SRPX, BTG3, ADK, TMEM51, ADM, PYGL, C17orf58, TRIM9, CSRP2, C1QTNF1, HEY1, TMSB4X, TIMP2, NEFL, ANXA7, SH3BGRL3, TIMP1, TEAD2, FGFBP1, SERPINB1, CHST7, ANXA1, ANXA2, GSTO1, IFNGR2, MYOF, ANXA5, DENND2A, EMP1, EMP3, DUSP6, SERPINB5, GCHFR, BACE2, TUBB2B, SLCO4A1, PXDN, PPARG, PLIN2, GRN, AHNAK, GSTP1, FSTL1, C11orf68, PAPSS2, SRPX2, CHN1, CLPTM1, HAS3, FLNC, CDKN2C, TPK1, CDKN2A, G0S2, DHRS3, DCBLD2, SULF2, FABP4, FABP5, DHRS9, ZYX, HSPA1A, TNC, PLOD2, CTGF, LAPTM4B, TUBB6, CDH3, CDH2, FCGRT, TUBB3, PLS3, IL13RA2, SLC16A3, TMSB10, PHLDA1, DGAT1, IGFBP4, FST, UBE2E1, RAB31, RAB34, MRTO4, FSCN1, S100A6, PKP1, S100A4, SQSTM1, FTL, RHOBTB3, KLC1, ITPRIPL2, PGK1, IGFBP6, PI3, PLTP, TFPI2, DRAP1, DAB2, TFCP2, P4HA1, TBXAS1, SPIRE1, OSBPL1A, NES, FHL1, FHL2, SGCE, SCRN1, MDK, EFHD2, PGM5, PGLS, UPP1, PGM1, ITGA3, LIX1L, ADORA2B, LCP1, MET, TNFRSF21, PPP1R15A, PCOLCE2, GLRX, ABLIM3, WBP5, CCL2, ZNF503, IGF2BP2, SH2B3, CTHRC1, SMOX, SEMA4B, AHNAK2, MYO1B, DBNDD2, MARCKS, POLR3A, EIF3G, OCIAD2, SEC24C, PFKP, ERRFI1, IFITM2, ANTXR2, HERC5, LGALS3, LGALS1, PEG10, FTH1, MAGEA1, TNS3, ST6GAL1, HLA-A, TM4SF1, HLA-E, BOLA2, TAGLN2, TAGLN3, TRIB1, RIN2, NBL1, CCNO, EPHA2, PFN2, SDC4, LPAR1, GSPT2, TSPAN5, DECR1, RPL39L, GNG11, PTK2, CIDEC, SPATS2L, CD24, ITM2C, NRP1, SLC35F2, SH3KBP1, MYL6B, PTPRK, AKAP12, UCHL1, RAC2, CD33, CTSC, CTSB, HTATIP2, ZCCHC24, TNFRSF12A, NME4, NOV, TNNT1, ARHGEF3, MAPRE3, FJX1, CD44, CSTB, FOXC1, GNAI2, SERTAD2, HSF1, HSD17B2, LMNA, CD55, MICB, CA12, CD70, CCDC109B, VEGFB, VEGFC, BOP1, KRT19, CD68, FOXA2, DFNA5, PRDM8, ETS1, ACTG1, TMEM145, RGS2, PLAU, ARHGDIB, COTL1, SIRPA, CD97, AMY1A, ANXA11, OCRL, ANXA10, KRT75, SCNN1G, PROCR, COL4A1, EHBP1L1, CHMP4A, VCL, RAI14, FBN2, HTRA1, PLD6, ADRB2, CACNA1H, IRAK2, TFAP2A, CRTAP, SMAD3, CAV2, TGFB3, CAV1, HMGA1, HMGA2, MSN, DNAJC12, ISG15, COL5A2, PHF19, DYRK2, TNFAIP6, CNRIP1, CXCL1, CXCL2, PITPNC1, HAPLN1, CA9, ABCC3, PLAUR, NAV2, ETV4, ETV5, TGFBR2, ALDH1A3, ACOX2, ALDH1A2, COL6A1, RBMS1, CUEDC1, FAM84B, AGPAT9, ASAP1, TFPI, AGPAT2, PDLIM1, C10orf11, CAMK2G, PDLIM7, SNCA, HES4, CMTM7, CMTM8, CMTM3, CARM1, ASB9, PMP22, GALM, VIM, CPVL".split(", ")));
		gen.setGenes(genes);
		gen.write("/home/babur/Documents/GeneSetEnrichment/Goutam/pathway");
		k.print();
	}
}
