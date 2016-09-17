package org.panda.misc.analyses;

import org.panda.misc.PrepareChiBENetworkDiffHighlightFile;
import org.panda.misc.PrepareChiBESIFNodeColorFormatFile;
import org.panda.resource.PCPathway;
import org.panda.utility.ArrayUtil;
import org.panda.utility.ValToColor;

import java.awt.*;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.*;

/**
 * @author Ozgun Babur
 */
public class CeciHumanAmnionProteomicsAnalysis
{
	public static final String BASE = "/home/babur/Documents/Analyses/Ceci/";
	public static void main(String[] args) throws IOException
	{
//		doEnrichment();
		prepareColorFile();
	}

	private static void doEnrichment() throws IOException
	{
		PCPathway.get().writeEnrichmentResults(genes, 3, 300, BASE + "HumanProtEnrichmentResults.txt");
	}

	private static void prepareColorFile() throws IOException
	{
		Map<String, Double> valMap = readValues();
		PrepareChiBESIFNodeColorFormatFile.prepare(valMap,
			new ValToColor(new double[]{-2, 0, 2}, new Color[]{new Color(100, 100, 255), Color.WHITE, new Color(255, 100, 100)}),
			BASE + "comparison.format");
	}

	private static Map<String, Double> readValues() throws FileNotFoundException
	{
		Scanner sc = new Scanner(new File("/home/babur/Documents/Analyses/Ceci/BRAC-223_hAm Nml Oligo_TMT_18fracs_20160706_OT_protein_rollup 071116.txt"));
		String[] header = sc.nextLine().split("\t");

		int geneInd = ArrayUtil.indexOf(header, "Gene");
		int valInd = ArrayUtil.indexOf(header, "FC");

		Map<String, Double> valMap = new HashMap<>();

		while (sc.hasNextLine())
		{
			String[] token = sc.nextLine().split("\t");
			String gene = token[geneInd].trim();
			if (gene.equals("DOY2")) gene = "DOPEY2";
			valMap.put(gene, Double.parseDouble(token[valInd]));
		}
		sc.close();
		return valMap;
	}

	private static Set<String> genes = new HashSet<>(Arrays.asList(
		("MFAP3L \n" +
		"MGST1 \n" +
		"SLPI \n" +
		"ZCCHC3 \n" +
		"DMBT1 \n" +
		"HMOX1 \n" +
		"NCCRP1 \n" +
		"ZG16B \n" +
		"MAOA \n" +
		"PLIN2 \n" +
		"ZBED6 \n" +
		"AOAH \n" +
		"SERPINB12 \n" +
		"ANP \n" +
		"AKAP12 \n" +
		"TAP2 \n" +
		"FZD1 \n" +
		"CKB \n" +
		"PAPPA2 \n" +
		"TAGLN \n" +
		"DOPEY2 \n" +
		"IGBP1 \n" +
		"MMRN1 \n" +
		"ASXL3 \n" +
		"COMMD1 \n" +
		"SLCO4A1 \n" +
		"SORD \n" +
		"ARHGEF16 \n" +
		"NPL \n" +
		"ADARB1 \n" +
		"MB \n" +
		"IGHD \n" +
		"TMCC2 \n" +
		"IGHA2 \n" +
		"PZP \n" +
		"TNC \n" +
		"AAK1 \n" +
		"TH \n" +
		"FAM129A \n" +
		"IFFO1 \n" +
		"FRS2 \n" +
		"CCZ1B \n" +
		"HBZ \n" +
		"ARFGAP2").replaceAll(" ", "").split("\n")));
}
