using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text;
using System.Text.RegularExpressions;
using System.Threading.Tasks;

namespace bioinf
{
	class Program
	{
		static string ViterbyToFindMostProbabilityPathBa10c(string str,
			char[] states,
			char[] letters,
			Dictionary<char, Dictionary<char, double>> emission,
			Dictionary<char, Dictionary<char, double>> transition)
		{
			var n = str.Length;

			var ptr = new Dictionary<char, char[]>();
			foreach (var state in states)
				foreach (var letter in letters)
					ptr[state] = new char[n];

			var matr = new Dictionary<char, List<double>>();
			foreach (var state in states)
			{
				matr[state] = new List<double>(n);
				for (var i = 0; i < n; i++)
					matr[state].Add(0);
				matr[state][0] = emission[state][str[0]] * 1 / states.Length;
			}

			double value;
			char key;
			for (var i = 1; i < n; i++)
			{
				foreach (var k in states)
				{
					value = 0;
					key = 'q';
					foreach (var state in states)
					{
						if (matr[state][i - 1] * transition[state][k] > value)
						{
							value = matr[state][i - 1] * transition[state][k];
							key = state;
						}
					}

					ptr[k][i] = key;
					matr[k][i] = emission[k][str[i]] * value;
				}

			}
			var answer = "";
			value = 0;
			var lastKey = 'q';
			foreach (var state in states)
			{
				if (matr[state][n - 1] > value)
				{
					value = matr[state][n - 1];
					lastKey = state;
				}
			}

			answer += lastKey;
			for (var i = n - 2; i >= 0; i--)
			{
				answer += ptr[lastKey][i + 1];
				lastKey = ptr[lastKey][i + 1];
			}
			answer = new string(answer.Reverse().ToArray());
			return answer;
		}

		static double ProbabilityOfStringByForwardAlgorithmBa10d(string str,
			char[] states,
			char[] letters,
			Dictionary<char, Dictionary<char, double>> emission,
			Dictionary<char, Dictionary<char, double>> transition)

		{
			var n = str.Length;
			var ptr = new Dictionary<char, char[]>();
			foreach (var state in states)
				foreach (var letter in letters)
					ptr[state] = new char[n];

			var matr = new Dictionary<char, List<double>>();
			foreach (var state in states)
			{
				matr[state] = new List<double>(n);
				for (var i = 0; i < n; i++)
					matr[state].Add(0);
				matr[state][0] = emission[state][str[0]] * 1 / states.Length;
			}

			double value;
			for (var i = 1; i < n; i++)
			{
				foreach (var k in states)
				{
					value = 0;
					foreach (var state in states)
						value += matr[state][i - 1] * transition[state][k];

					matr[k][i] = emission[k][str[i]] * value;
				}

			}
			value = 0;
			foreach (var state in states)
				value += matr[state][n - 1];

			return value;
		}

		static Dictionary<string, Dictionary<string, double>> CreateEmptyTransitionMatrix(int strLen)
		{
			var transition = new Dictionary<string, Dictionary<string, double>>();
			transition.Add("S", new Dictionary<string, double>());
			transition.Add("I0", new Dictionary<string, double>());
			for (var i = 1; i <= strLen; i++)
			{
				transition.Add($"M{i}", new Dictionary<string, double>());
				transition.Add($"D{i}", new Dictionary<string, double>());
				transition.Add($"I{i}", new Dictionary<string, double>());
			}
			transition.Add("E", new Dictionary<string, double>());

			foreach (var value in transition.Values)
			{
				value.Add("S", 0);
				value.Add("I0", 0);
				for (var i = 1; i <= strLen; i++)
				{
					value.Add($"M{i}", 0);
					value.Add($"D{i}", 0);
					value.Add($"I{i}", 0);
				}
				value.Add("E", 0);
			}
			return transition;
		}

		static List<List<int>> SplitInsertionColumnNumbers(List<int> notSplitted)
		{
			if (notSplitted.Count == 0)
			{
				return new List<List<int>>();
			}
			var splitted = new List<List<int>>();
			splitted.Add(new List<int>());
			splitted[0].Add(notSplitted[0]);
			foreach (var i in notSplitted.Skip(1))
			{
				if (splitted.Last().Last() + 1 == i)
					splitted.Last().Add(i);
				else
				{
					splitted.Add(new List<int>());
					splitted.Last().Add(i);
				}
			}
			return splitted;
		}

		static bool ContainsOnly(string str, char c)
		{
			for (var i = 0; i < str.Length; i++)
			{
				if (str[i] != c)
					return false;
			}
			return true;
		}

		static List<int> IndexesOfNotInsert(List<List<char>> insertionColumns)
		{
			var n = insertionColumns[0].Count;
			var transpose = new List<string>();
			for (var i = 0; i < n; i++)
			{
				transpose.Add(new string(insertionColumns.Select(s => s[i]).ToArray()));
			}

			var indexes = new List<int>();
			for (var i = 0; i < n; i++)
			{
				if (ContainsOnly(transpose[i], '-'))
					indexes.Add(i);

			}
			return indexes;
		}

		static List<int> IndexesOfNotMatch(List<char> matchColumn)
		{
			var indexes = new List<int>();
			for (var i = 0; i < matchColumn.Count; i++)
			{
				if (matchColumn[i] == '-')
					indexes.Add(i);

			}
			return indexes;
		}

		static int CountOfToInsert(List<List<char>> insertionColumns)
		{
			var n = insertionColumns[0].Count;
			var transpose = new List<string>();
			for (var i = 0; i < n; i++)
			{
				transpose.Add(new string(insertionColumns.Select(s => s[i]).ToArray()));
			}

			var count = 0;
			for (var i = 0; i < n; i++)
			{
				count += Math.Max(transpose[i].Count(x => x != '-') - 1, 0);

			}
			return count;
		}

		static List<TransitionDataColumn> CreateInnerViewOfAlignedStrings(List<string> alignedStrings, double threshold)
		{
			var n = alignedStrings.Count;
			var strLen = alignedStrings[0].Length;

			var transposeStrings = new List<string>();
			for (var i = 0; i < strLen; i++)
			{
				transposeStrings.Add(new string(alignedStrings.Select(s => s[i]).ToArray()));
			}

			var data = new List<TransitionDataColumn>();
			data.Add(new TransitionDataColumn());
			for (var i = 0; i < strLen; i++)
			{
				if (transposeStrings[i].Count(x => x == '-') / Convert.ToDouble(n) < threshold)
				{
					data.Add(new TransitionDataColumn());
					data.Last().Match = transposeStrings[i].ToList();
				}
				else
				{
					data.Last().Insertion.Add(transposeStrings[i].ToList());
				}
			}
			return data;
		}

		static Dictionary<string, Dictionary<string, double>> CreateTransitionMatrix(List<string> alignedStrings, double threshold)
		{
			var n = alignedStrings.Count;
			var strLen = alignedStrings[0].Length;

			var data = CreateInnerViewOfAlignedStrings(alignedStrings, threshold);

			var fullRange = new List<int>();
			for (var i = 0; i < n; i++)
				fullRange.Add(i);

			var transition = CreateEmptyTransitionMatrix(data.Count - 1);

			List<int> indexesOfNotMatch;
			IEnumerable<int> indexesOfMatch;
			List<int> indexesOfNotInsert;
			IEnumerable<int> indexesOfInsert;
			int goToInsert;
			List<int> indexesOfNotMatchI1;
			IEnumerable<int> indexesOfMatchI1;
			for (var i = 1; i < data.Count-1; i++)
			{
				indexesOfNotMatch = IndexesOfNotMatch(data[i].Match);
				indexesOfNotMatchI1 = IndexesOfNotMatch(data[i + 1].Match);
				indexesOfMatch = fullRange.Except(indexesOfNotMatch);
				indexesOfMatchI1 = fullRange.Except(indexesOfNotMatchI1);

				if (data[i].Insertion.Count != 0)
				{
					indexesOfNotInsert = IndexesOfNotInsert(data[i].Insertion);
					indexesOfInsert = fullRange.Except(indexesOfNotInsert);

					transition[$"M{i}"][$"I{i}"] = Convert.ToDouble(indexesOfMatch.Except(indexesOfNotInsert).Count())/
					                               indexesOfMatch.Count();
					transition[$"M{i}"][$"D{i + 1}"] = Convert.ToDouble(indexesOfMatch.Except(indexesOfInsert).Except(indexesOfMatchI1).Count())/
													   indexesOfMatch.Count();
					transition[$"M{i}"][$"M{i + 1}"] = Convert.ToDouble(indexesOfMatch.Except(indexesOfInsert).Except(indexesOfNotMatchI1).Count())/
													   indexesOfMatch.Count();

					goToInsert = CountOfToInsert(data[i].Insertion);
					var goToMatch = indexesOfInsert.Except(indexesOfNotMatchI1).Count();
					var goToDelete = indexesOfInsert.Except(indexesOfMatchI1).Count();
					if (goToInsert + goToMatch + goToDelete == 0)
						continue;
					transition[$"I{i}"][$"I{i}"] = Convert.ToDouble(goToInsert)/
					                               (goToInsert + goToMatch + goToDelete);
					transition[$"I{i}"][$"D{i + 1}"] = Convert.ToDouble(goToDelete)/
					                                   (goToInsert + goToMatch + goToDelete);
					transition[$"I{i}"][$"M{i + 1}"] = Convert.ToDouble(goToMatch)/
					                                   (goToInsert + goToMatch + goToDelete);

					goToInsert = indexesOfNotMatch.Except(indexesOfNotInsert).Count();
					goToDelete = indexesOfNotMatch.Except(indexesOfInsert).Except(indexesOfMatchI1).Count();
					goToMatch = indexesOfNotMatch.Except(indexesOfInsert).Except(indexesOfNotMatchI1).Count();
					if (goToInsert + goToMatch + goToDelete == 0)
						continue;
					transition[$"D{i}"][$"I{i}"] = Convert.ToDouble(goToInsert)/
					                               (goToInsert + goToMatch + goToDelete);
					transition[$"D{i}"][$"D{i + 1}"] = Convert.ToDouble(goToDelete)/
					                                   (goToInsert + goToMatch + goToDelete);
					transition[$"D{i}"][$"M{i + 1}"] = Convert.ToDouble(goToMatch)/
					                                   (goToInsert + goToMatch + goToDelete);
				}
				else
				{
					var goToDelete = indexesOfMatch.Except(indexesOfMatchI1).Count();
					var goToMatch = indexesOfMatch.Except(indexesOfNotMatchI1).Count();
					if (goToDelete + goToMatch == 0)
						continue;
					transition[$"M{i}"][$"D{i + 1}"] = Convert.ToDouble(goToDelete)/(goToDelete + goToMatch);
					transition[$"M{i}"][$"M{i + 1}"] = Convert.ToDouble(goToMatch)/ (goToDelete + goToMatch);

					goToDelete = indexesOfNotMatch.Except(indexesOfMatchI1).Count();
					goToMatch = indexesOfNotMatch.Except(indexesOfNotMatchI1).Count();
					if (goToDelete + goToMatch == 0)
						continue;
					transition[$"D{i}"][$"D{i + 1}"] = Convert.ToDouble(goToDelete)/(goToDelete + goToMatch);
					transition[$"D{i}"][$"M{i + 1}"] = Convert.ToDouble(goToMatch)/(goToDelete + goToMatch);
				}
			}

			var j = 0;
			indexesOfNotMatchI1 = IndexesOfNotMatch(data[1].Match);
			indexesOfMatchI1 = fullRange.Except(indexesOfNotMatchI1);
			if (data[j].Insertion.Count != 0)
			{
				indexesOfNotInsert = IndexesOfNotInsert(data[j].Insertion);
				indexesOfInsert = fullRange.Except(indexesOfNotInsert);
				goToInsert = indexesOfNotInsert.Count();
				var goToDelete = indexesOfNotInsert.Except(indexesOfMatchI1).Count();
				var goToMatch = indexesOfNotInsert.Except(indexesOfNotMatchI1).Count();
				if (goToInsert + goToMatch + goToDelete != 0)
				{
					transition[$"E"][$"I{0}"] = Convert.ToDouble(goToInsert)/
					                            (goToInsert + goToMatch + goToDelete);
					transition[$"E"][$"D{1}"] = Convert.ToDouble(goToDelete)/
					                            (goToInsert + goToMatch + goToDelete);
					transition[$"E"][$"M{1}"] = Convert.ToDouble(goToMatch)/
					                            (goToInsert + goToMatch + goToDelete);
				}

				goToInsert = CountOfToInsert(data[j].Insertion);
				goToDelete = indexesOfInsert.Except(indexesOfMatchI1).Count();
				goToMatch = indexesOfInsert.Except(indexesOfNotMatchI1).Count();
				if (goToInsert + goToMatch + goToDelete != 0)
				{
					transition[$"I{0}"][$"I{0}"] = Convert.ToDouble(goToInsert)/
					                               (goToInsert + goToMatch + goToDelete);
					transition[$"I{0}"][$"D{1}"] = Convert.ToDouble(goToDelete)/
					                               (goToInsert + goToMatch + goToDelete);
					transition[$"I{0}"][$"M{1}"] = Convert.ToDouble(goToMatch)/
					                               (goToInsert + goToMatch + goToDelete);
				}
			}
			else
			{
				var goToDelete = indexesOfNotMatchI1.Count();
				var goToMatch = indexesOfMatchI1.Count();
				transition[$"E"][$"D{1}"] = Convert.ToDouble(goToDelete)/n;
				transition[$"E"][$"M{1}"] = Convert.ToDouble(goToMatch)/n;
			}

			j = data.Count - 1;
			indexesOfNotMatch = IndexesOfNotMatch(data[j].Match);
			indexesOfMatch = fullRange.Except(indexesOfNotMatch);
			indexesOfNotInsert = IndexesOfNotInsert(data[j].Insertion);
			indexesOfInsert = fullRange.Except(indexesOfNotInsert);
			if (data[j].Insertion.Count != 0)
			{
				goToInsert = indexesOfMatch.Except(indexesOfNotInsert).Count();
				var goToE = indexesOfMatch.Except(indexesOfInsert).Count();
				if (goToInsert + goToE != 0)
				{
					transition[$"M{j}"][$"I{j}"] = Convert.ToDouble(goToInsert) / (goToInsert + goToE);
					transition[$"M{j}"][$"E"] = Convert.ToDouble(goToE) / (goToInsert + goToE);
				}

				goToInsert = indexesOfNotMatch.Except(indexesOfNotInsert).Count();
				goToE = indexesOfNotMatch.Except(indexesOfInsert).Count();
				if (goToInsert + goToE != 0)
				{
					transition[$"D{j}"][$"I{j}"] = Convert.ToDouble(goToInsert) / (goToInsert + goToE);
					transition[$"D{j}"][$"E"] = Convert.ToDouble(goToE) / (goToInsert + goToE);
				}

				goToInsert = CountOfToInsert(data[j].Insertion);
				goToE = indexesOfInsert.Count();
				if (goToInsert + goToE != 0)
				{
					transition[$"I{j}"][$"I{j}"] = Convert.ToDouble(goToInsert) / (goToInsert + goToE);
					transition[$"I{j}"][$"E"] = Convert.ToDouble(goToE) / (goToInsert + goToE);
				}
			}
			else
			{
				transition[$"M{j}"][$"E"] = Convert.ToDouble(indexesOfMatch.Count()) / n;
				transition[$"D{j}"][$"E"] = Convert.ToDouble(indexesOfNotMatch.Count()) / n;
			}

			return transition;
		}

		static Dictionary<string, Dictionary<char, double>> CreateEmptyEmissionMatrix(int strLen, List<char> alphabet)
		{
			var emission = new Dictionary<string, Dictionary<char, double>>();
			emission.Add("S", new Dictionary<char, double>());
			emission.Add("I0", new Dictionary<char, double>());
			for (var i = 1; i <= strLen; i++)
			{
				emission.Add($"M{i}", new Dictionary<char, double>());
				emission.Add($"D{i}", new Dictionary<char, double>());
				emission.Add($"I{i}", new Dictionary<char, double>());
			}
			emission.Add("E", new Dictionary<char, double>());

			foreach (var value in emission.Values)
			{
				foreach (var c in alphabet)
				{
					value.Add(c, 0);
				}
				
			}
			return emission;
		}

		static Dictionary<string, Dictionary<char, double>> CreateEmissionMatrix(
			List<string> alignedStrings, 
			double threshold,
			List<char> alphabet)
		{
			var n = alignedStrings.Count;
			var strLen = alignedStrings[0].Length;

			var data = CreateInnerViewOfAlignedStrings(alignedStrings, threshold);

			var emission = CreateEmptyEmissionMatrix(data.Count - 1, alphabet);

			return emission;

		}

		static void Main(string[] args)
		{
			var strings = "DCDABACED.DCCA--CA-.DCDAB-CA-.BCDA---A-.BC-ABE-AE";
			var threshold = 0.252;
            var matrix = CreateTransitionMatrix(strings.Split('.').ToList(), threshold);
			using (StreamWriter writetext = new StreamWriter("write.txt"))
			{
				writetext.Write($" \t");
				foreach (var key in matrix.Keys)
				{
					writetext.Write($"{key}\t");
				}
				writetext.WriteLine();

				foreach (var str in matrix)
				{
					writetext.Write($" {str.Key}\t");
					foreach (var val in str.Value)
					{
						writetext.Write($"{val.Value}\t");
					}
					writetext.WriteLine();
				}
			}

			var emission = CreateEmissionMatrix(strings.Split('.').ToList(), threshold, new List<char>() {'A', 'B', 'C', 'D', 'E'});

		}
	}
}
