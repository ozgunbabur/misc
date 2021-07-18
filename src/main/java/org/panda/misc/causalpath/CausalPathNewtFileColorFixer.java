package org.panda.misc.causalpath;

import org.panda.utility.FileUtil;

import java.io.BufferedWriter;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.List;
import java.util.stream.Collectors;

public class CausalPathNewtFileColorFixer
{
	static String[] files = ("Molm14-parental-vs-early.nwt\n" +
		"Molm14-parental-vs-late.nwt\n" +
		"Molm14-early-vs-late.nwt\n" +
		"MV411-early-vs-parental.nwt\n" +
		"MV411-late-vs-early.nwt\n" +
		"MV411-late-vs-parental.nwt.nwt").split("\n");

	public static void main(String[] args) throws IOException
	{
		for (String file : files)
		{
			fix("/Users/ozgun/Downloads/" + file);
		}
	}


	public static void fix(String file) throws IOException
	{
		List<String> content = Files.lines(Paths.get(file)).collect(Collectors.toList());

		BufferedWriter writer = Files.newBufferedWriter(Paths.get(file));

		content.forEach(l ->
		{
			l = l.replaceAll("#009700", "#009600");
			l = l.replaceAll("#62392d", "#009600");
			l = l.replaceAll("#8c1005", "#960000");
			FileUtil.writeln(l, writer);
		});

		writer.close();
	}
}
