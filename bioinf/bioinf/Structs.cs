using System.Collections.Generic;

namespace bioinf
{
	public class TransitionDataColumn
	{
		public List<char> Match;
		public List<List<char>> Insertion;

		public TransitionDataColumn()
		{
			Match = new List<char>();
			Insertion = new List<List<char>>();
		}
	}
}