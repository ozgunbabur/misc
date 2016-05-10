package org.panda.misc.rppa.parser;

import org.cbio.causality.network.PhosphoSitePlus;
import org.cbio.causality.rppa.RPPAData;
import org.cbio.causality.rppa.RPPAFileReader;
import org.cbio.causality.rppa.RPPANetworkMapper;
import org.panda.utility.CollectionUtil;
import org.panda.utility.ValToColor;

import java.awt.*;
import java.io.*;
import java.util.*;
import java.util.List;

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
		int valInd1 = ACHN;
		int valInd2 = ACHN;

		Set<RPPAData> datas = new HashSet<>();
		datas.addAll(convertRPPA(base + "PP_data_Table1_Thomas_Q177800_pY1000_1PC_FDR_121913.csv", valInd1));
		datas.addAll(convertRPPA(base + "PP_data_Table2_Thomas_Q177800_Basophillic_1PC_FDR_121913.csv", valInd2));
		datas.addAll(readRPPAData(valInd1 == ACHN));

		RPPAData mtor = new RPPAData("MTOR-inactivated", null, Arrays.asList("MTOR"), null);
		mtor.makeActivityNode(false);
		datas.add(mtor);

		RPPAData.ChangeAdapter chDet = new RPPAData.ChangeAdapter(){};

		double threshold = 1.5;
		chDet.setThreshold(threshold);
		for (RPPAData data : datas) if (!data.isActivity()) data.setChDet(chDet);

//		prepareChiBEFormatFile(datas, base + "achn.format");
//		if (true) return;

		// Fill-in missing effect from PhosphoSitePlus
		PhosphoSitePlus.fillInMissingEffect(datas, 10);

		RPPANetworkMapper.writeGraph(datas, threshold, base + "ACHN-compatible.sif",
			RPPANetworkMapper.GraphType.COMPATIBLE, null);

//		writeRPPA(data, outFile);
	}

	private static void compareSets(Set<RPPAData>... sets)
	{
		Set<String>[] ids = new Set[sets.length];

		for (int i = 0; i < ids.length; i++)
		{
			ids[i] = new HashSet<>();
			for (RPPAData data : sets[i])
			{
				ids[i].add(data.id);
			}
		}
		CollectionUtil.printVennSets(ids);
	}

	private static Set<RPPAData> convertRPPA(String inFile, int valindex) throws FileNotFoundException
	{
		Set<RPPAData> set = new HashSet<>();
		Scanner sc = new Scanner(new File(inFile));
		while (sc.hasNextLine())
		{
			RPPAData data = readLine(sc.nextLine(), valindex);
			if (data != null) set.add(data);
		}
		return set;
	}

	private static void writeRPPA(Set<RPPAData> set, String outFile)
	{
		List<RPPAData> list = new ArrayList<>(set);
		Collections.sort(list, (o1, o2) -> o1.id.compareTo(o2.id));
		RPPAData.write(list, outFile);
	}

	private static void prepareChiBEFormatFile(Set<RPPAData> datas, String filename) throws IOException
	{
		Map<String, Double> map = new HashMap<>();
		for (RPPAData data : datas)
		{
			if (data.isActivity()) continue;

			for (String gene : data.genes)
			{
				double v = Math.abs(data.getChangeValue());
				if (!map.containsKey(gene) || map.get(gene) < v) map.put(gene, v);
			}
		}

		ValToColor v2c = new ValToColor(new double[]{1, 5}, new Color[]{Color.WHITE, Color.GREEN});

		BufferedWriter writer = new BufferedWriter(new FileWriter(filename));
		writer.write("node\tall-nodes\tcolor\t255 255 255");
		writer.write("\nnode\tall-nodes\tbordercolor\t0 0 0");

		for (String gene : map.keySet())
		{
			double v = map.get(gene);
			Color c = v2c.getColor(v);
			writer.write("\nnode\t" + gene + "\tcolor\t" + c.getRed() + " " + c.getGreen() + " " + c.getBlue());
			if (v >= 1.5) writer.write("\nnode\t" + gene + "\tborderwidth\t2");
			writer.write("\nnode\t" + gene + "\ttooltip\t" + v);
		}
		
		writer.close();
	}
	
	private static RPPAData readLine(String line, int valIndex)
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

		double[][] vals = new double[0][];

		if (valIndex >= 0)
		{
			double v = Double.parseDouble(token[valIndex]);
			vals = new double[1][1];
			vals[0][0] = v;
		}
		else if (valIndex == TABLE1_CONTROLS)
		{
			if (token[36].isEmpty() || token[38].isEmpty()) return null;

			double vA = Double.parseDouble(token[36]);
			double vS = Double.parseDouble(token[38]);

			double foldCh = Math.max(vA, vS) / Math.min(vA, vS);
			if (vS < vA) foldCh *= -1;
			vals = new double[1][1];
			vals[0][0] = foldCh;
		}
		else if (valIndex == TABLE2_CONTROLS)
		{
			if (token[34].isEmpty() || token[36].isEmpty()) return null;

			double vA = Double.parseDouble(token[34]);
			double vS = Double.parseDouble(token[36]);

			double foldCh = Math.max(vA, vS) / Math.min(vA, vS);
			if (vS < vA) foldCh *= -1;
			vals = new double[1][1];
			vals[0][0] = foldCh;
		}

		return new RPPAData(generateID(geneList, siteMap), vals, geneList, siteMap);
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

	private static Set<RPPAData> readRPPAData(boolean achn) throws FileNotFoundException
	{
		List<RPPAData> datas = RPPAFileReader.readAnnotation(base + "abdata-chibe.txt", "ID1", "Symbols", "Sites", "Effect");

		Map<String, RPPAData> map = new HashMap<>();
		for (RPPAData data : datas)
		{
			String id = data.id.substring(0, data.id.indexOf("_GBL"));
			map.put(id, data);
		}

		Set<RPPAData> set = new HashSet<>();
		Scanner sc = new Scanner(new File(base + "RPPA_data_Thomas_lab_SR_18Feb2016.csv"));
		sc.nextLine();

		while (sc.hasNextLine())
		{
			String[] token = sc.nextLine().split("\t");
			String id = token[0].replaceAll("-", ".");
			id = id.substring(0, id.indexOf("_GBL"));

			RPPAData data = map.get(id);
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

			data.vals = new double[][]{{fold}};
			set.add(data);
		}

		sc.close();
		return set;
	}


}
