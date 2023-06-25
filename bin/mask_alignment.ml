open! Base

let version = "1.0.0"

let usage_msg =
  {eof|usage: mask_alignment <aln.fa> <max_gap_percent> > <output.fa>

example: mask_alignment silly.aln.fa 95 > silly.aln.masked.fa
         ^^^ That would make fasta file that has any column with >= 95% gaps removed.

Note: To get rid of columns with all gaps, pass in 100.  
|eof}

let help_msg = [%string "mask_alignment version %{version}\n\n%{usage_msg}"]

let abort ?(exit_code = 1) msg =
  let () = Stdio.eprintf "%s\n" msg in
  Caml.exit exit_code

let parse_mask_ratio s =
  match Float.of_string s with
  | exception exn ->
      let msg = Exn.to_string exn in
      abort [%string "ERROR -- mask_ratio can't be coerced to Float: %{msg}"]
  | x ->
      if Float.(x <= 0. || x > 100.) then
        abort [%string "ERROR -- max_gap_percent must be in the range (0, 100]"]
      else x /. 100.

let gap_ratio num_gaps num_seqs = Int.to_float num_gaps /. Int.to_float num_seqs

(* Read in the records. Check that all alignment lengths are good. *)
let get_records fname =
  let open Bio_io.Fasta in
  (* let aln_len, num_seqs, records = *)
  let aln_len, num_seqs, records =
    In_channel.with_file_foldi_records fname ~init:(0, 0, [])
      ~f:(fun i (aln_len, num_seqs, records) record ->
        let open Record in
        let seq_len = String.length @@ seq @@ record in
        let expected_aln_len = if Int.(i = 0) then seq_len else aln_len in
        if seq_len <> expected_aln_len then
          abort
            [%string
              "ERROR -- %{id record} should be %{expected_aln_len#Int} bases \
               but was %{seq_len#Int} bases!"] ;
        (expected_aln_len, num_seqs + 1, record :: records) )
  in
  (aln_len, num_seqs, Array.of_list_rev records)

let a = Char.to_int 'a'
let z = Char.to_int 'z'
let a_cap = Char.to_int 'A'
let z_cap = Char.to_int 'Z'

let is_gap c =
  let c' = Char.to_int c in
  not ((a <= c' && c' <= z) || (a_cap <= c' && c' <= z_cap))

let get_good_columns ~aln_len ~num_seqs records max_gap_ratio =
  let keep_these = ref [] in
  let () =
    for column_i = 0 to aln_len - 1 do
      (* Number of gaps in this column. *)
      let num_gaps = ref 0 in
      let () =
        for seq_i = 0 to num_seqs - 1 do
          let record = records.(seq_i) in
          let seq = Bio_io.Fasta.Record.seq record in
          let char = String.get seq column_i in
          if is_gap char then num_gaps := !num_gaps + 1
        done
      in
      if Float.(gap_ratio !num_gaps num_seqs < max_gap_ratio) then
        keep_these := column_i :: !keep_these
    done
  in
  Array.of_list_rev !keep_these

(* Parse input args. *)
let fname, max_gap_ratio =
  match Sys.get_argv () with
  | [|_; fname; max_gap_ratio|] ->
      (fname, parse_mask_ratio max_gap_ratio)
  | _ ->
      abort help_msg

let aln_len, num_seqs, records = get_records fname
let good_columns = get_good_columns ~aln_len ~num_seqs records max_gap_ratio
let masked_seq_length = Array.length good_columns

let print_masked_alignment () =
  let open Bio_io.Fasta in
  Array.iter records ~f:(fun record ->
      let seq = record |> Record.seq in
      (* Reusing this buffer saves a bit of time, but not worth it. Simpler to
         just keep it right here. *)
      let buf = Bytes.create masked_seq_length in
      (* The loop is a smidge faster, but less nice... *)
      for i = 0 to masked_seq_length - 1 do
        let char = String.get seq good_columns.(i) in
        Bytes.set buf i char
      done ;
      let new_seq = Bytes.to_string buf in
      let new_record = Record.with_seq new_seq record in
      Stdio.print_endline @@ Record.to_string new_record )

let () = print_masked_alignment ()
