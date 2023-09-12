open! Base

let abort ?(exit_code = 1) msg =
  let () = Stdio.eprintf "%s\n" msg in
  Caml.exit exit_code

module Alignment : sig
  (** - [num_cols] is the number of columns in the alignment (aka alignment
        length) *)
  type t = private {num_cols: int; records: Bio_io.Fasta.Record.t array}

  val num_seqs : t -> int
  (** [num_seqs t] returns the number of sequences in the alignment. *)

  val read_alignment_file : string -> t
  (** [read_alignment_file file_name] reads the alignment specified by
      [file_name], returning a valid [t].

      It assures that the [num_cols] and [num_seqs] are correct. Thus they can
      be used with [unsafe_get] functions downstream, if that consumer takes a
      [t] value.

      Raises an exception unless all the records have the same length (number
      alignment columns). *)
end = struct
  type t = {num_cols: int; records: Bio_io.Fasta.Record.t array}

  let num_seqs t = Array.length t.records

  let read_alignment_file fname =
    let open Bio_io.Fasta in
    let num_cols, records =
      In_channel.with_file_foldi_records fname ~init:(0, [])
        ~f:(fun i (num_cols, records) record ->
          let open Record in
          let seq_len = String.length @@ seq @@ record in
          let expected_num_cols = if Int.(i = 0) then seq_len else num_cols in
          if seq_len <> expected_num_cols then
            abort
              [%string
                "ERROR -- %{id record} should be %{expected_num_cols#Int} \
                 bases but was %{seq_len#Int} bases!"] ;
          (expected_num_cols, record :: records) )
    in
    {num_cols; records= Array.of_list_rev records}
end

let version = "1.0.0"

let usage_msg =
  {eof|usage: mask_alignment <aln.fa> <max_gap_percent> > <output.fa>

example: mask_alignment silly.aln.fa 95 > silly.aln.masked.fa
         ^^^ That would make fasta file that has any column with >= 95% gaps removed.

Note: To get rid of columns with all gaps, pass in 100.  
|eof}

let help_msg = [%string "mask_alignment version %{version}\n\n%{usage_msg}"]

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

let a = Char.to_int 'a'
let z = Char.to_int 'z'
let a_cap = Char.to_int 'A'
let z_cap = Char.to_int 'Z'

let is_gap c =
  let c' = Char.to_int c in
  not ((a <= c' && c' <= z) || (a_cap <= c' && c' <= z_cap))

let get_good_columns (Alignment.{num_cols; records} as alignment) max_gap_ratio
    =
  let keep_these = ref [] in
  let num_seqs = Alignment.num_seqs alignment in
  let () =
    for column_i = 0 to num_cols - 1 do
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

let alignment = Alignment.read_alignment_file fname
let good_columns = get_good_columns alignment max_gap_ratio
let masked_seq_length = Array.length good_columns

let print_masked_alignment () =
  let open Bio_io.Fasta in
  alignment.records
  |> Array.iter ~f:(fun record ->
         let seq = record |> Record.seq in
         (* Reusing this buffer saves a bit of time, but not worth it. Simpler
            to just keep it right here. *)
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
