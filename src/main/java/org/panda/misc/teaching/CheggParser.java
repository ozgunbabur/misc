package org.panda.misc.teaching;

import org.panda.utility.CollectionUtil;
import org.panda.utility.FileUtil;

import java.io.BufferedWriter;
import java.io.IOException;
import java.util.List;
import java.util.Set;
import java.util.stream.Collectors;

public class CheggParser
{
	public static void main(String[] args) throws IOException
	{
		String dir = "/Users/ozgun/Documents/Teaching/CS220/grades-spring-2021/final-cheating/Chegg";
		writeHTMLs(dir + "/AskerDetail.csv", dir);
	}

	public static void writeHTMLs(String inFile, String outDir) throws IOException
	{
		Set<String> users = getUsers(inFile);

		for (String user : users)
		{
			BufferedWriter writer = FileUtil.newBufferedWriter(outDir + "/" + user + ".html");

			for (int i = 1; i <= 7; i++)
			{
				String question = "Question " + i;
				writer.write("<p>" + question + "</p>\n");
				writer.write(getData(inFile, user, question) + "\n");
			}

			writer.close();
		}

	}

	public static Set<String> getUsers(String file)
	{
		return FileUtil.getTermsInTabDelimitedColumn(file, 4, 1).stream().distinct().map(s -> s.substring(0, s.indexOf("@"))).collect(Collectors.toSet());
	}

	public static String getData(String file, String user, String questionInd)
	{
		List<String> list = FileUtil.linesTabbed(file).skip(1).filter(t -> t[4].contains(user) && t[8].contains(questionInd)).map(t -> outOfQuotes(t[9])).collect(Collectors.toList());
		return CollectionUtil.merge(list, "\n");
	}

	public static String outOfQuotes(String s)
	{
		if (s.startsWith("\"")) s = s.substring(1);
		if (s.endsWith("\"")) s = s.substring(0, s.length() - 1);
		return s.replaceAll("\"\"", "\"");
	}
}
