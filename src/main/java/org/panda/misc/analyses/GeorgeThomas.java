package org.panda.misc.analyses;

import org.biopax.paxtools.pattern.miner.SIFEnum;
import org.panda.causalpath.resource.ProteomicsFileReader;
import org.panda.causalpath.run.CausalityAnalysisSingleMethodInterface;
import org.panda.resource.PhosphoSitePlus;
import org.panda.resource.network.PathwayCommons;
import org.panda.resource.network.SignedPC;
import org.panda.resource.signednetwork.SignedType;
import org.panda.resource.tcga.ProteomicsFileRow;
import org.panda.utility.ArrayUtil;
import org.panda.utility.CollectionUtil;
import org.panda.utility.graph.DirectedGraph;
import org.panda.utility.graph.Graph;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.*;
import java.util.stream.Collectors;

/**
 * Created by babur on 3/2/16.
 */
public class GeorgeThomas
{
	public static final int ACHN = 2;
	public static final int SN12C = 3;
	public static final int TABLE1_CONTROLS = -5;
	public static final int TABLE2_CONTROLS = -6;
	public static final Set<String> MTOR_RELATED = new HashSet<>(Arrays.asList("AKT1S1", "AKT2", "AKT3", "BAD",
		"EIF4EBP1", "IRS2", "MAF1", "PDCD4", "RPS6KB1", "TBC1D4", "MTOR", "DEPDC6", "DAP1", "AKT1", "SGK", "PIK3CA",
		"ULK1", "LPIN1", "ATG13", "NFAT3", "STAT3", "CCT", "EEF2K", "EIF4B", "GSK3A", "GSK3B", "MAD1", "RICTOR",
		"POLDIP3", "BRAF", "CASP9", "CDKN1A", "CDKN1B", "CHEK1", "CHUK", "CRAF", "FOXO1", "FOXO3", "FOXO4", "MAP3K5",
		"NOS3", "PPARGC1A", "RAF1", "SKP2", "TSC2", "WNK1",
		"EIF4EBP2", "MDM2", "NDRG1", "NDRG2", "RPS6", "RAPTOR"));

	static String base = "/home/babur/Documents/RPPA/GeorgeThomas_phosphoproteomics/";

	public static void main(String[] args) throws IOException
	{
		doCausalityAnalysis();
//		reportTFEnrichment(readRNAseq());
	}

	public static void doCausalityAnalysis() throws IOException
	{
		int valInd1 = ACHN;
		int valInd2 = ACHN;

		Set<ProteomicsFileRow> rppas = new HashSet<>();
		rppas.addAll(convertRPPA(base + "PP_data_Table1_Thomas_Q177800_pY1000_1PC_FDR_121913.csv", valInd1));
		rppas.addAll(convertRPPA(base + "PP_data_Table2_Thomas_Q177800_Basophillic_1PC_FDR_121913.csv", valInd2));
		rppas.addAll(readProteomicsFileRow(valInd1 == ACHN));

		// Fill-in missing effect from PhosphoSitePlus
		PhosphoSitePlus.get().fillInMissingEffect(rppas, 10);

		ProteomicsFileRow mtor = new ProteomicsFileRow("MTOR-inactivated", null, Arrays.asList("MTOR"), null);
		mtor.makeActivityNode(true);
		rppas.add(mtor);

		CausalityAnalysisSingleMethodInterface.generateCausalityGraph(rppas, 1.5, "compatible", true, 0, true, 1,
			base + "ACHN-compatible-temp2");
	}

	private static void compareSets(Set<ProteomicsFileRow>... sets)
	{
		Set<String>[] ids = new Set[sets.length];

		for (int i = 0; i < ids.length; i++)
		{
			ids[i] = new HashSet<>();
			for (ProteomicsFileRow data : sets[i])
			{
				ids[i].add(data.id);
			}
		}
		CollectionUtil.printVennSets(ids);
	}

	private static Set<ProteomicsFileRow> convertRPPA(String inFile, int valindex) throws FileNotFoundException
	{
		Set<ProteomicsFileRow> set = new HashSet<>();
		Scanner sc = new Scanner(new File(inFile));
		while (sc.hasNextLine())
		{
			ProteomicsFileRow data = readLine(sc.nextLine(), valindex);
			if (data != null) set.add(data);
		}
		return set;
	}

	private static ProteomicsFileRow readLine(String line, int valIndex)
	{
		String[] token = line.split("\t");

		if (token.length < 11) return null;
		if (token[10].isEmpty()) return null;
		if (valIndex >= 0 && token[valIndex].isEmpty()) return null;
		if (valIndex >= 0 && token[valIndex].equals("N.D.")) return null;

		// Phosphorylation exists but site do not match with known sequence.
		if (token[20].contains("*") && token[12].isEmpty()) return null;

		String[] genes = token[10].split("; ");
		String[] prots = token[11].split("; ");
		String[] sites = token[12].replaceAll("%", "").split("; ");

		List<String> aa = getMarkedAminoAcids(token[20]);
		List<String> geneList = new ArrayList<>();
		Map<String, List<String>> siteMap = new HashMap<>();

		assert genes.length == prots.length && genes.length == sites.length;

		for (int i = 0; i < genes.length; i++)
		{
			if (prots[i].contains(" iso")) continue;

			geneList.add(genes[i]);

			String[] pos = sites[i].split(", ");

			if (pos.length != aa.size() && aa.size() != 1)
			{
				System.out.println(line);
				return null;
			}

			List<String> siteList = new ArrayList<>();
			for (int j = 0; j < pos.length; j++)
			{
				String a = aa.size() == 1 ? aa.get(0) : aa.get(j);
				siteList.add(a +  pos[j]);
			}
			siteMap.put(genes[i], siteList);
		}

		if (geneList.isEmpty()) return null;

		double[] vals = new double[0];

		if (valIndex >= 0)
		{
			double v = Double.parseDouble(token[valIndex]);
			vals = new double[]{v};
		}
		else if (valIndex == TABLE1_CONTROLS)
		{
			if (token[36].isEmpty() || token[38].isEmpty()) return null;

			double vA = Double.parseDouble(token[36]);
			double vS = Double.parseDouble(token[38]);

			double foldCh = Math.max(vA, vS) / Math.min(vA, vS);
			if (vS < vA) foldCh *= -1;
			vals = new double[]{foldCh};
		}
		else if (valIndex == TABLE2_CONTROLS)
		{
			if (token[34].isEmpty() || token[36].isEmpty()) return null;

			double vA = Double.parseDouble(token[34]);
			double vS = Double.parseDouble(token[36]);

			double foldCh = Math.max(vA, vS) / Math.min(vA, vS);
			if (vS < vA) foldCh *= -1;
			vals = new double[]{foldCh};
		}

		return new ProteomicsFileRow(generateID(geneList, siteMap), vals, geneList, siteMap);
	}

	private static String generateID(List<String> genes, Map<String, List<String>> sites)
	{
		StringBuilder sb = new StringBuilder();

		boolean first = true;
		for (String gene : genes)
		{
			if (first)
			{
				sb.append(gene);
				first = false;
			}
			else sb.append("--").append(gene);

			for (String s : sites.get(gene))
			{
				sb.append("-").append(s);
			}
		}
		return sb.toString();
	}

	private static List<String> getMarkedAminoAcids(String s)
	{
		List<String> list = new ArrayList<>();

		while (s.contains("*"))
		{
			int i = s.indexOf("*");
			String aa = s.substring(i - 1, i);
			list.add(aa);
			s = s.substring(i + 1);
		}
		return list;
	}

	private static Set<ProteomicsFileRow> readProteomicsFileRow(boolean achn) throws FileNotFoundException
	{
		List<ProteomicsFileRow> datas = ProteomicsFileReader.readAnnotation(base + "abdata-chibe.txt", "ID1", "Symbols", "Sites", "Effect");

		Map<String, ProteomicsFileRow> map = new HashMap<>();
		for (ProteomicsFileRow data : datas)
		{
			String id = data.id.substring(0, data.id.indexOf("_GBL"));
			map.put(id, data);
		}

		Set<ProteomicsFileRow> set = new HashSet<>();
		Scanner sc = new Scanner(new File(base + "RPPA_data_Thomas_lab_SR_18Feb2016.csv"));
		sc.nextLine();

		while (sc.hasNextLine())
		{
			String[] token = sc.nextLine().split("\t");
			String id = token[0].replaceAll("-", ".");
			id = id.substring(0, id.indexOf("_GBL"));

			ProteomicsFileRow data = map.get(id);
			if (data == null)
			{
				System.out.println("No match found = " + token[0]);
				continue;
			}

			double v0 = Double.parseDouble(token[achn ? 1 : 3]);
			double v1 = Double.parseDouble(token[achn ? 2 : 4]);

			v0 = Math.pow(2, v0);
			v1 = Math.pow(2, v1);

			double fold = Math.max(v0, v1) / Math.min(v0, v1);
			if (v1 < v0) fold = -fold;

			data.vals = new double[]{fold};
			set.add(data);
		}

		sc.close();
		return set;
	}

	private static Map<String, Double> readRNAseq() throws IOException
	{
		String filename = base + "RNAseq_ACHN_CYT_updatedTable3_March172014-1.txt";
		String[] header = Files.lines(Paths.get(filename)).skip(3).findFirst().get().split("\t");
		int geneIndex = ArrayUtil.indexOf(header, "Gene Symbol");
		int fcIndex = ArrayUtil.indexOf(header, "FC_Cyt24_Veh24");

		return Files.lines(Paths.get(filename)).skip(4).map(l -> l.split("\t"))
			.collect(Collectors.toMap(t -> t[geneIndex], t -> Double.parseDouble(t[fcIndex]),
				(v1, v2) -> Math.abs(v1) > Math.abs(v2) ? v1 : v2));
	}

	private static void reportTFEnrichment(Map<String, Double> rnaseq)
	{
		double thr = 2;
		long posCnt = rnaseq.values().stream().filter(v -> Math.abs(v) >= thr).filter(v -> v > 0).count();
		long negCnt = rnaseq.values().stream().filter(v -> Math.abs(v) >= thr).filter(v -> v < 0).count();

		System.out.println("posCnt = " + posCnt);
		System.out.println("negCnt = " + negCnt);

		DirectedGraph graph = (DirectedGraph) PathwayCommons.get().getGraph(SIFEnum.CONTROLS_EXPRESSION_OF);
		Map<String, Double> select = rnaseq.keySet().stream().filter(k -> Math.abs(rnaseq.get(k)) >= thr)
			.collect(Collectors.toMap(k -> k, rnaseq::get));

		List<String> enriched = graph.getEnrichedGenes(
			select.keySet(), rnaseq.keySet(), 0.1, DirectedGraph.NeighborType.DOWNSTREAM, 1, 5);

		System.out.println("enriched = " + enriched);

		DirectedGraph upGraph = SignedPC.get().getGraph(SignedType.UPREGULATES_EXPRESSION);
		DirectedGraph dwGraph = SignedPC.get().getGraph(SignedType.DOWNREGULATES_EXPRESSION);

		for (String gene : enriched)
		{
			System.out.println("\ngene = " + gene + "\t" + graph.getDownstream(gene).size() + "\t" +
				CollectionUtil.countOverlap(graph.getDownstream(gene), rnaseq.keySet()) + "\t" +
				CollectionUtil.countOverlap(graph.getDownstream(gene), select.keySet()));

			System.out.println("Upregulated\t" + CollectionUtil.countOverlap(upGraph.getDownstream(gene), rnaseq.keySet()) + "\t" + CollectionUtil.countOverlap(upGraph.getDownstream(gene), select.keySet()));
			for (String target : upGraph.getDownstream(gene))
			{
				if (select.containsKey(target)) System.out.println(target + "\t" + rnaseq.get(target));
			}
			System.out.println("Downregulated	" + CollectionUtil.countOverlap(dwGraph.getDownstream(gene), rnaseq.keySet()) + "\t" + CollectionUtil.countOverlap(dwGraph.getDownstream(gene), select.keySet()));
			for (String target : dwGraph.getDownstream(gene))
			{
				if (select.containsKey(target)) System.out.println(target + "\t" + rnaseq.get(target));
			}
		}
	}
}
