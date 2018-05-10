package org.panda.misc.analyses;

import org.panda.utility.FileUtil;

import java.io.BufferedWriter;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;

/**
 * @author Ozgun Babur
 */
public class Vyse
{
	public static final String DIR = "/home/babur/Documents/Analyses/Vyse-drug-resistance/";

	public static void main(String[] args) throws IOException
	{
		convert(DIR + "PazR-vs-parental-A204.txt", DIR + "PazR/data.txt");
		convert(DIR + "DazR-vs-parental-A204.txt", DIR + "DazR/data.txt");
	}

	public static void convert(String inFile, String outFile) throws IOException
	{
		BufferedWriter writer = Files.newBufferedWriter(Paths.get(outFile));
		writer.write("ID\tSymbols\tSites\tEffect\tValue");
		Files.lines(Paths.get(inFile)).skip(2).map(l -> l.split("\t"))
			.filter(t -> t.length > 11 && (t[10].equals("+") || t[11].equals("+"))).forEach(t ->
		{
			String id = t[0];
			String gene = t[3];
			String site = t[16] + t[17];

			double v = -Math.log(Double.valueOf(t[8]));
			if (t[9].startsWith("-")) v *= -1;

			FileUtil.lnwrite(id + "\t" + gene + "\t" + site + "\t\t" + v, writer);
		});
		writer.close();
	}
}
