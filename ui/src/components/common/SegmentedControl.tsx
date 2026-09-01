import { useId } from 'react';

export interface SegmentOption<T extends string> {
  value: T;
  label: string;
  hint?: string;
}

/** Radiogroup-styled segmented control. Keyboard: arrows move + select. */
export function SegmentedControl<T extends string>({
  legend,
  name,
  options,
  value,
  onChange,
  describedById,
}: {
  legend: string;
  name: string;
  options: SegmentOption<T>[];
  value: T;
  onChange: (value: T) => void;
  describedById?: string;
}) {
  const groupId = useId();
  return (
    <fieldset
      className="segmented"
      role="radiogroup"
      aria-labelledby={`${groupId}-legend`}
      aria-describedby={describedById}
    >
      <legend id={`${groupId}-legend`} className="field-label">
        {legend}
      </legend>
      <div className="segmented__track">
        {options.map((opt) => {
          const id = `${groupId}-${opt.value}`;
          const checked = opt.value === value;
          return (
            <label
              key={opt.value}
              htmlFor={id}
              className={`segmented__seg${checked ? ' is-selected' : ''}`}
            >
              <input
                type="radio"
                id={id}
                name={name}
                value={opt.value}
                checked={checked}
                onChange={() => onChange(opt.value)}
                className="sr-only"
              />
              <span className="segmented__label">{opt.label}</span>
              {opt.hint ? <span className="segmented__hint">{opt.hint}</span> : null}
            </label>
          );
        })}
      </div>
    </fieldset>
  );
}
