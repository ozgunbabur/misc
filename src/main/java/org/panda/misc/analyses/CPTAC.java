package org.panda.misc.analyses;

import org.panda.resource.UniProtSequence;
import org.panda.utility.ArrayUtil;
import org.panda.utility.CollectionUtil;
import org.panda.utility.FileUtil;

import java.io.BufferedWriter;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.List;

public class CPTAC
{
	public static void main(String[] args) throws IOException
	{
		fixAFile("PTRC_baseline_combined_with_stoich_20181029_exp3.txt");
		fixAFile("PTRC_baseline_combined_with_stoich_20190318_exp2.txt");
		fixAFile("ptrc_ex8_4plex_combined_with_stoich_20190426_exp8.txt");
		fixAFile("ptrc_ex11_kurtz_GlobalAndPhospho_ColCoded_exp11.txt");
	}

	public static void fixAFile(String name) throws IOException
	{
		String dir = "/home/ozgun/Analyses/Druker/datasets/";
		fillInMissingSites(dir + name, dir + "fixed/" + name);
	}

	public static void fillInMissingSites(String inFile, String outFile) throws IOException
	{
		BufferedWriter writer = Files.newBufferedWriter(Paths.get(outFile));

		Files.lines(Paths.get(inFile)).forEach(l ->
		{
			if (l.startsWith("Gene\tEntry"))
			{
				FileUtil.writeln(l, writer);
			}
			else
			{
				String[] t = l.split("\t");

				if (!t[3].isEmpty() || t[4].equals("NA"))
				{
					FileUtil.writeln(l, writer);
				}
				else
				{
					String up = t[1];
					String pep = t[4];

					String wStars = pep.replaceAll("\\.", "").replaceAll("-", "");
					String plain = wStars.replaceAll("\\*", "");

					int start = UniProtSequence.get().getStartLocation(up, plain);

					if (start > 0)
					{
						List<String> sites = new ArrayList<>();
						String[] x = wStars.split("\\*");
						int base = start;
						for (int i = 0; i < x.length - 1; i++)
						{
							String s = x[i];
							String aa = s.substring(s.length() - 1);
							sites.add(aa + (base + s.length() - 1));

							base += s.length();
						}
						t[3] = CollectionUtil.merge(sites, "|");
						FileUtil.writeln(ArrayUtil.merge("\t", t), writer);
					}
				}
			}
		});


		writer.close();
	}
}
