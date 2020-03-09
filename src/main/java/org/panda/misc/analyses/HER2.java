package org.panda.misc.analyses;

import org.panda.misc.CausalPathSubnetwork;
import org.panda.utility.FileUtil;
import org.panda.utility.SIFFileUtil;
import org.panda.utility.StreamDirection;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.Arrays;
import java.util.HashSet;

import static org.panda.utility.StreamDirection.BOTHSTREAM;

public class HER2
{
	public static final String DIR = "/Users/ozgun/Documents/Analyses/CP-paper-runs/TCGA-RPPA/";

	public static void main(String[] args) throws IOException
	{
//		printAllEGFR();
		generateRecurrenceSubgraph();
	}

	private static void printAllEGFR() throws IOException
	{
		Files.list(Paths.get(DIR)).forEach(p ->
		{
			System.out.println("p = " + p);
			FileUtil.lines(p + "/causative.sif")
				.filter(l -> l.contains("EGFR") || l.contains("ERBB2") || l.contains("ERBB3"))
				.forEach(System.out::println);
		});
	}

	private static void generateRecurrenceSubgraph() throws IOException
	{
		SIFFileUtil.writeNeighborhood(DIR + "/recurrent/1.sif", new HashSet<>(Arrays.asList("EGFR", "ERBB2", "ERBB3")), DIR + "/recurrent/1-EGFR.sif", BOTHSTREAM);
	}
}
