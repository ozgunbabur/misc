package org.panda.misc.analyses;

import org.biopax.paxtools.pattern.miner.SIFEnum;
import org.panda.misc.PathwayEnrichmentSIFGenerator;
import org.panda.resource.HGNC;
import org.panda.resource.PCPathway;
import org.panda.utility.CollectionUtil;
import org.panda.utility.Kronometre;
import org.panda.utility.statistics.FDR;

import java.awt.*;
import java.io.*;
import java.util.*;
import java.util.List;
import java.util.stream.Collectors;

import static java.util.stream.Collectors.*;

/**
 * Created by babur on 2/29/16.
 */
public class Goutam
{
	static String DIR = "/home/babur/Documents/Analyses/Goutam/";

	public static void main(String[] args) throws IOException
	{
		Kronometre k = new Kronometre();

//		Set<String> smad = getBMPSMADPathwayGenes();
//		System.out.println("smad.size() = " + smad.size());
//		PCPathway.get().addCustomPathway("bmp.smad.pathway", "BMP-SMAD Signaling", smad, Collections.emptySet());

		Set<String> allGenes = getUnion(readAllGenes(0));
		Set<String> upGenes = getUnion(readAllGenes(1));
		Set<String> dwGenes = getUnion(readAllGenes(-1));
		Set<String> undecidedGenes = CollectionUtil.getIntersection(upGenes, dwGenes);

		Set<String> pathwayIDs = doEnrichment(1);
		Set<String> set = doEnrichment(-1);
		if (set != null) pathwayIDs.addAll(set);

//		if (true) return;

		PathwayEnrichmentSIFGenerator sg = new PathwayEnrichmentSIFGenerator();
		sg.setMolecules(allGenes);
		sg.setOwlFilename("/home/babur/Documents/PC/PathwayCommons.8.Detailed.BIOPAX.owl");
		sg.setBlacklistFilename("/home/babur/Documents/PC/blacklist.txt");
		sg.setPathwayIDs(pathwayIDs);
		String[][] groupTerms = {
			new String[]{"collagen", "extracellular"},
			new String[]{"integrin"},
			new String[]{"AP-1", "AP1"},
//			new String[]{"TGF-beta", "TGF_beta"},
//			new String[]{"syndecan"},
//			new String[]{"fibrin"},
//			new String[]{"RA biosynthesis", "retinoid", "visual"},
//			new String[]{"C-MYC"},
//			new String[]{"HIF-1"},
//			new String[]{"Scavenging"},
//			new String[]{"VEGF"},
//			new String[]{"FGF"},
//			new String[]{"miRNA"},
//			new String[]{"pyrimidine"},
//			new String[]{"Creatine"}
		};

		sg.setGroupTerms(groupTerms);
		sg.setTypes(SIFEnum.CONTROLS_STATE_CHANGE_OF, SIFEnum.CONTROLS_EXPRESSION_OF);

		Color upColor = new Color(255, 180, 200);
		Color dwColor = new Color(180, 200, 255);
		Color mixColor = new Color(255, 255, 220);
		allGenes.forEach(g -> sg.addNodeColor(g,
			undecidedGenes.contains(g) ? mixColor : upGenes.contains(g) ? upColor : dwColor));

		sg.write(DIR + "pathway");
		k.print();
	}

	private static Set<String> doEnrichment(int sign) throws IOException
	{
		System.out.println("sign = " + sign);
		List<Set<String>> sets = readAllGenes(sign);
		Set<String> genes = getIntersection(sets);
		showOverlaps(sets);
		return findEnrichedPathways(genes, sign);
	}

	private static void showOverlaps(List<Set<String>> bag) throws IOException
	{
		CollectionUtil.printNameMapping("Caso Cl 2", "MDV Cl 1", "Neo Sh#5", "Neo Sh# 3");
		CollectionUtil.printVennSets(bag.toArray(new Collection[bag.size()]));
	}

	private static Set<String> findEnrichedPathways(Set<String> genes, Integer sign) throws IOException
	{
		String filename = DIR + "results" + (sign == 0 ? "" : sign > 0 ? "-upregulated" : "-downregulated") + ".txt";

		int minMember = 3;
		int maxMember = 300;
		PCPathway.get().writeEnrichmentResults(genes, minMember, maxMember, filename);
		Map<String, Double>[] pvals = PCPathway.get().getEnrichmentPvals(genes, null, minMember, maxMember);
		Map<String, Double> qvals = FDR.getQVals(pvals[0], pvals[1]);
		OptionalDouble thr = pvals[0].keySet().stream().filter(id -> qvals.get(id) < 0.1).mapToDouble(pvals[0]::get).max();
		if (!thr.isPresent()) return null;

		return pvals[0].keySet().stream().filter(id -> pvals[0].get(id) <= thr.getAsDouble()).collect(toSet());
	}

	private static Set<String> getIntersection(List<Set<String>> bag)
	{
		Set<String> intersection = new HashSet<>();

		boolean first = true;

		for (Set<String> set : bag)
		{
			if (first)
			{
				intersection.addAll(set);
				first = false;
			}
			else
			{
				intersection.retainAll(set);
			}
		}
		return intersection;
	}

	private static Set<String> getUnion(List<Set<String>> bag)
	{
		return bag.stream().flatMap(Collection::stream).collect(Collectors.toSet());
	}

	private static List<Set<String>> readAllGenes(Integer sign) throws FileNotFoundException
	{
		List<Set<String>> bag = new ArrayList<>();
		bag.add(readGenes("/home/babur/Downloads/AnalysisList/Caso Cl 2 vs Contol copy.csv", sign));
		bag.add(readGenes("/home/babur/Downloads/AnalysisList/MDV Cl 1 vs Control copy.csv", sign));
		bag.add(readGenes("/home/babur/Downloads/AnalysisList/Neo Sh#5  vs Control copy.csv", sign));
		bag.add(readGenes("/home/babur/Downloads/AnalysisList/Neo Sh# 3 vs Control copy.csv", sign));
		return bag;
	}

	private static Set<String> readGenes(String file, int sign) throws FileNotFoundException
	{
		Set<String> set = new HashSet<>();
		Scanner sc = new Scanner(new File(file));
		sc.nextLine();
		while (sc.hasNextLine())
		{
			String[] token = sc.nextLine().split("\t");
			String gene = token[2];
			if (gene.isEmpty()) continue;
			String foldChange = token[5];

			if (sign == 0 || (sign < 0 && foldChange.startsWith("-")) || sign > 0 && !foldChange.startsWith("-"))
				set.add(gene);
		}
		return set;
	}

	private static Set<String> getBMPSMADPathwayGenes()
	{
		String s = "Foxi3\n" +
			"Cacna1g\n" +
			"Bhlhe41\n" +
			"Nav2\n" +
			"Fcgbp\n" +
			"Foxo6\n" +
			"Fcgr3\n" +
			"Calml4\n" +
			"C1qc\n" +
			"Edn2\n" +
			"Slc22a4\n" +
			"Hsd11b2\n" +
			"C1qc\n" +
			"Sned1\n" +
			"Lgals9\n" +
			"Cav1\n" +
			"Tgfbr3\n" +
			"Igfbp4\n" +
			"Nxph3\n" +
			"Qrfp\n" +
			"Cyth4\n" +
			"Fbln2\n" +
			"Adamts17\n" +
			"Ccdc164\n" +
			"Angptl2\n" +
			"Filip1l\n" +
			"Fcgbp\n" +
			"Scube1\n" +
			"Ngfr\n" +
			"Lgr6\n" +
			"Pld4\n" +
			"Tril\n" +
			"Nrxn2\n" +
			"Lrig1\n" +
			"Chdh\n" +
			"Lmo2\n" +
			"Edar\n" +
			"Inhbb\n" +
			"Pde4b\n" +
			"Sema5a\n" +
			"Aif1l\n" +
			"Pdgfb\n" +
			"Parvg\n" +
			"Scmh1\n" +
			"Sox9\n" +
			"Mcf2l\n" +
			"Rasgrp2\n" +
			"Fam189a2\n" +
			"Nbl1\n" +
			"Arc\n" +
			"Gm973\n" +
			"Ggta1\n" +
			"Nav2\n" +
			"Lmo1\n" +
			"Ptprv\n" +
			"Krt19\n" +
			"Trp73\n" +
			"Necab3\n" +
			"Lhx2\n" +
			"Sema3g\n" +
			"Col16a1\n" +
			"Arhgef16\n" +
			"Epas1\n" +
			"Fam20c\n" +
			"Spry4\n" +
			"Wnt9a\n" +
			"Sorcs2\n" +
			"Lrp4\n" +
			"Flrt1\n" +
			"Bcl11b\n" +
			"Cplx2\n" +
			"Meox1\n" +
			"Tspan18\n" +
			"Rab11fip4\n" +
			"Pdlim4\n" +
			"Irx4\n" +
			"D0H4S114\n" +
			"Crabp2\n" +
			"Phlda1\n" +
			"Barx2\n" +
			"Mfap3l\n" +
			"Sorcs2\n" +
			"Il34\n" +
			"Rgs11\n" +
			"Eya2\n" +
			"Shisa2\n" +
			"Cyp46a1\n" +
			"Fry\n" +
			"Ypel4\n" +
			"Ptch2\n" +
			"Dnmt3a\n" +
			"Fcgbp\n" +
			"D0H4S114\n" +
			"Slc26a6\n" +
			"Sec14l2\n" +
			"Lbh\n" +
			"2310003H01Rik\n" +
			"Ece1\n" +
			"Ephb4\n" +
			"Lynx1\n" +
			"Card14\n" +
			"Pmepa1\n" +
			"Has3\n" +
			"9930012K11Rik\n" +
			"Anxa2\n" +
			"Sema3b\n" +
			"Rapgef3\n" +
			"Iffo2\n" +
			"Mesp2\n" +
			"Card10\n" +
			"Mef2a\n" +
			"Foxq1\n" +
			"Gpr153\n" +
			"Nrp2\n" +
			"Slco3a1\n" +
			"Dlx2\n" +
			"Dbn1\n" +
			"Ptges\n" +
			"Tpm2\n" +
			"Ttc39a\n" +
			"Spint1\n" +
			"Prrt1\n" +
			"P2ry2\n" +
			"Cited4\n" +
			"Dmpk\n" +
			"Msx2\n" +
			"Pla2g2e\n" +
			"Kcnh1\n" +
			"Hes2\n" +
			"Engase\n" +
			"Ovol1\n" +
			"Id2\n" +
			"Mt4\n" +
			"Dlx4\n" +
			"Hspb8\n" +
			"Id3\n" +
			"A730090H04Rik\n" +
			"Cldn4\n" +
			"Ccdc64b\n" +
			"Krt35\n" +
			"Krt82\n" +
			"Hephl1\n" +
			"Acpl2\n" +
			"S100a3\n" +
			"Vsig8\n" +
			"Sult2b1\n" +
			"5430421N21Rik\n" +
			"Gng13\n" +
			"Ablim2";

		System.out.println("s.split(\"\\n\").length = " + s.split("\n").length);
		return Arrays.stream(s.split("\n")).map(g -> HGNC.get().getSymbol(s)).filter(Objects::nonNull).collect(toSet());
	}
}
